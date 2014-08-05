/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
  *
 * Copyright (c) 2014 University of Basel
 * All rights reserved.
 */


#include "itkBoundingBox.h"
#include "itkStandardMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkPenalizingMeanSquaresPointSetToImageMetric.h"
#include "itkLBFGSOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkCommand.h"
#include "itkMesh.h"
#include "itkTransformMeshFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkPosteriorModelBuilder.h"
#include "itkPointsLocator.h"
#include "itkPointSetToImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "utils.h"


//
// config values for the optimizer
//
struct Config {
    // how much margin will be added around the bounding box of the reference shape
    const static double imageMargin;

    // exit criterion for the optimizer
    static const double gradientTolerance;

    // the number of iterations for the optimization
    static const unsigned numberOfIterations;
};
const double Config::imageMargin = 50;
const double Config::gradientTolerance = 1e-12;
const unsigned Config::numberOfIterations = 100;


//
// the types that we need later on
//
const unsigned Dimensions = 3;
typedef itk::PointSet<float, Dimensions  > PointSetType;
typedef itk::Mesh<float, Dimensions  > MeshType;
typedef itk::Point<double, 3> PointType;
typedef itk::Image<float, Dimensions> DistanceImageType;
typedef itk::StandardMeshRepresenter<float, Dimensions> RepresenterType;
typedef itk::PenalizingMeanSquaresPointSetToImageMetric<PointSetType, DistanceImageType> MetricType;
typedef itk::StatisticalShapeModelTransform<MeshType, double, Dimensions> StatisticalModelTransformType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::PointSetToImageRegistrationMethod<PointSetType,DistanceImageType > RegistrationFilterType;
typedef itk::LBFGSOptimizer OptimizerType;
typedef itk::Image< unsigned char, Dimensions > BinaryImageType;
typedef itk::DanielssonDistanceMapImageFilter<BinaryImageType, DistanceImageType> DistanceFilterType;
typedef itk::TransformMeshFilter<MeshType, MeshType, StatisticalModelTransformType> TransformMeshFilterType;
typedef itk::PointSetToImageFilter<PointSetType,BinaryImageType> PointsToImageFilterType;
typedef itk::LinearInterpolateImageFunction<DistanceImageType, double> InterpolatorType;

typedef itk::PosteriorModelBuilder<MeshType> PosteriorModelBuilderType;
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
typedef itk::PointsLocator< MeshType::PointsContainer > PointsLocatorType;
#else
typedef itk::PointsLocator<int, 3, double, MeshType::PointsContainer > PointsLocatorType;
#endif


 /*
 * Compute a new model, which is restricted to go through the proints specified in targetLandmarks.
 * The variance argument specifies the assumed landmark inaccuracy (assuming the landmark inaccuracy is modeled as a zero-mean gaussian).
 */
StatisticalModelType::Pointer
computePosteriorModel( const StatisticalModelType* statisticalModel,
                       const  std::vector<PointType >& modelLandmarks,
                      const  std::vector<PointType >& targetLandmarks,
                      double variance)
{

        StatisticalModelType::PointValueListType constraints;

        // We need to make sure the the points in fixed landmarks are real vertex points of the model reference.
        MeshType::Pointer reference = statisticalModel->GetRepresenter()->GetReference();
        PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
        ptLocator->SetPoints(reference->GetPoints());
        ptLocator->Initialize();

        assert(modelLandmarks.size() == targetLandmarks.size());
        for (unsigned i = 0; i < targetLandmarks.size(); i++) {

            int closestPointId = ptLocator->FindClosestPoint(modelLandmarks[i]);
            PointType refPoint = (*reference->GetPoints())[closestPointId];

            // compensate for the rigid transformation that was applied to the model
            StatisticalModelType::PointValuePairType pointValue(refPoint ,targetLandmarks[i]);
            constraints.push_back(pointValue);

        }

        PosteriorModelBuilderType::Pointer PosteriorModelBuilder = PosteriorModelBuilderType::New();
        StatisticalModelType::Pointer PosteriorModel = PosteriorModelBuilder->BuildNewModelFromModel(statisticalModel,constraints, variance, false);

        return PosteriorModel;
}


/*
 * Compute a distance image from the given mesh. The distance image will be the size of the bounding box of the point set,
 * plus the given margin.
 */
DistanceImageType::Pointer distanceImageFromMesh(MeshType* mesh, double margin) {

    // Compute a bounding box around the reference shape
    typedef  itk::BoundingBox<int, 3, float, MeshType::PointsContainer> BoundingBoxType;
    BoundingBoxType::Pointer bb = BoundingBoxType::New();
    bb->SetPoints(mesh->GetPoints());
    bb->ComputeBoundingBox();

    // Compute a binary image from the point set, which is as large as the bounding box plus a margin.
    PointsToImageFilterType::Pointer pointsToImageFilter = PointsToImageFilterType::New();
    pointsToImageFilter->SetInput( mesh );
    BinaryImageType::SpacingType spacing; spacing.Fill( 1.0 );
    BinaryImageType::PointType origin = bb->GetMinimum();
    BinaryImageType::SpacingType diff = bb->GetMaximum() - bb->GetMinimum();
    BinaryImageType::SizeType size;
    for (unsigned i =0; i < 3; i++) {
        origin[i] -= margin; // a five cm margin
        size[i] = diff[i] + margin;
    }


    pointsToImageFilter->SetSpacing( spacing );
    pointsToImageFilter->SetOrigin( origin );
    pointsToImageFilter->SetSize( size);
    pointsToImageFilter->Update();

    // compute a distance map to the points in the pointset
    BinaryImageType::Pointer binaryImage = pointsToImageFilter->GetOutput();
    DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
    distanceFilter->SetInput( binaryImage );
    distanceFilter->Update();
    DistanceImageType::Pointer distanceImage = distanceFilter->GetOutput();
    return distanceImage;
}



/*
 * Project every point of mesh to the closest point on the target and returns the
 * resulting projected (mesh).
 */
MeshType::Pointer projectOnTargetMesh(MeshType* mesh, MeshType* targetMesh) {

    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(targetMesh->GetPoints());
    ptLocator->Initialize();

    MeshType::Pointer projectedMesh = Utils::cloneMesh<MeshType>(mesh);
    for (unsigned i = 0; i < mesh->GetNumberOfPoints(); i++) {
        MeshType::PointType pt = mesh->GetPoint(i);
        int closestPointId = ptLocator->FindClosestPoint(pt);
        MeshType::PointType closestPt = targetMesh->GetPoint(closestPointId);
        projectedMesh->SetPoint(i, closestPt);
    }
    return projectedMesh;

}

/*
 * This class is used to track the progress of the optimization
 * (its method Execute is called in each iteration of the optimization)
*/
class IterationStatusObserver : public itk::Command
{
public:
  typedef  IterationStatusObserver   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;

  itkNewMacro( Self );

  typedef itk::LBFGSOptimizer    OptimizerType;
  //typedef itk::RegularStepGradientDescentBaseOptimizer OptimizerType;

  typedef const OptimizerType                     *OptimizerPointer;


   void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
    {
      return;
    }

    std::cout << "Iteration: " << ++m_iter_no ;
   std::cout << "; Value: " << optimizer->GetCachedValue();
    std::cout << "; Current Parameters: " << optimizer->GetCachedCurrentPosition() << std::endl;
  }


protected:
  IterationStatusObserver():
     m_iter_no(0)     {};

  virtual ~IterationStatusObserver(){};

private:
  int m_iter_no;

};

/*
 * Perform a registration (fitting) of the Gaussian Process model to the target mesh. Write the resulting mesh into the file specified
 * by outputmesh. Optionally, the user can specify matching landmark points on the shape. These will be matched during the fitting.
 * The lmvariance argument specifies the assumed landmark inaccuracy (assuming the landmark inaccuracy is modeled as a zero-mean gaussian).
 */
int main(int argc, char* argv[]) {


    if (argc < 6) {
        std::cout << "usage " << argv[0] << " modelname targetmesh regWeight fittedmesh projectedMesh [fixedLandmarks movingLandmarks lmVariance] " << std::endl;
        exit(-1);
    }

    char* modelName = argv[1];
    char* targetName = argv[2];
    double regWeight = atof(argv[3]);
    char* fittedMeshName = argv[4];
    char* projectedMeshName = argv[5];

    char* fixedLandmarksName = 0;
    char* movingLandmarksName = 0;
    double lmVariance = 0;

    if (argc == 9) {
        fixedLandmarksName = argv[6];
        movingLandmarksName = argv[7];
        lmVariance = atof(argv[8]);
    }



    // load the Surface to which we will fit
    MeshType::Pointer targetMesh = Utils::readMesh<MeshType>(targetName);

    // load the model
    StatisticalModelType::Pointer model = StatisticalModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    model->Load(representer, modelName);

    StatisticalModelType::Pointer constrainedModel = 0;

    // in case that landmarks were defined, we constrain the model on the landmarks
    if (fixedLandmarksName != 0 && movingLandmarksName != 0) {

        std::cout << "constraining model on landmarks" << std::endl;
        // read the landmarks
        std::vector<PointType> fixedLandmarks = Utils::readLandmarks(fixedLandmarksName);
        std::vector<PointType> movingLandmarks = Utils::readLandmarks(movingLandmarksName);

        // compute a new model, which is constrained to match the landmark points
        constrainedModel = computePosteriorModel(model, fixedLandmarks, movingLandmarks, lmVariance);
    }
    else {
        constrainedModel = model; // use the original model
    }


    // create a shape model transform
    StatisticalModelTransformType::Pointer statModelTransform = StatisticalModelTransformType::New();
    statModelTransform->SetStatisticalModel(constrainedModel);
    statModelTransform->SetIdentity();


    // The actual fitting will be done to a distance image representation of the mesh.
    DistanceImageType::Pointer distanceImage = distanceImageFromMesh(targetMesh, Config::imageMargin);


    // set up the optimizer
    OptimizerType::Pointer optimizer = OptimizerType::New();
    const unsigned long numberOfIterations = Config::numberOfIterations;
    const double gradientTolerance = Config::gradientTolerance; // convergence criterion

    optimizer->SetMaximumNumberOfFunctionEvaluations(numberOfIterations );
    optimizer->SetGradientConvergenceTolerance(gradientTolerance );
    optimizer->MinimizeOn();

    // set up the observer to keep track of the progress
    typedef  IterationStatusObserver ObserverType;
    ObserverType::Pointer observer = ObserverType::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );

    // set up the metric and interpolators
    MetricType::Pointer metric = MetricType::New();
    metric->SetRegularizationParameter(regWeight);
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // we create the fixedPointSet of the registration from the reference mesh of our model.
    // As we are fitting to the 0 level set of a distance image, we set the value of the pointdata to 0.
    PointSetType::Pointer fixedPointSet = PointSetType::New();
    fixedPointSet->SetPoints(model->GetRepresenter()->GetReference()->GetPoints());
    PointSetType::PointDataContainer::Pointer points = PointSetType::PointDataContainer::New();
    points->Reserve(model->GetRepresenter()->GetReference()->GetNumberOfPoints());
    for (PointSetType::PointDataContainer::Iterator it = points->Begin(); it != points->End(); ++it) {
        it->Value() = 0;
    }
    fixedPointSet->SetPointData(points);


    // connect all the components
    RegistrationFilterType::Pointer registration = RegistrationFilterType::New();
    registration->SetInitialTransformParameters(statModelTransform->GetParameters());
    registration->SetMetric(metric);
    registration->SetInterpolator(interpolator);
    registration->SetOptimizer(   optimizer);
    registration->SetTransform(   statModelTransform );

    registration->SetFixedPointSet(  fixedPointSet);
    registration->SetMovingImage(distanceImage);

    try {
        std::cout << "starting model fitting" << std::endl;
        registration->Update();

    } catch ( itk::ExceptionObject& o ) {
        std::cout << "caught exception " << o << std::endl;
    }


    // Transform the mesh according to the optimized transformation.
    TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
    transformMeshFilter->SetInput(model->GetRepresenter()->GetReference());
    transformMeshFilter->SetTransform(statModelTransform);
    transformMeshFilter->Update();


    // We might be interested in both the fitting result (which is an approximation to the target
    // or an actual projection of the fitted mesh onto the target. We compute and  write both.

    MeshType::Pointer fittedMesh = transformMeshFilter->GetOutput();
    MeshType::Pointer projectedMesh = projectOnTargetMesh(fittedMesh, targetMesh);

    Utils::writeMesh<MeshType>(fittedMesh, fittedMeshName);
    Utils::writeMesh<MeshType>(projectedMesh, projectedMeshName);

}



