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
#include "itkMeanSquaresPointSetToImageMetric.h"
#include "itkLBFGSOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkCommand.h"
#include "itkMesh.h"
#include "itkTransformMeshFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkPosteriorModelBuilder.h"
#include "itkPointsLocator.h"
#include "itkPointSetToImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

const unsigned Dimensions = 3;
typedef itk::PointSet<float, Dimensions  > PointSetType;
typedef itk::Mesh<float, Dimensions  > MeshType;
typedef itk::Point<double, 3> PointType;
typedef itk::Image<float, Dimensions> DistanceImageType;
typedef itk::StandardMeshRepresenter<float, Dimensions> RepresenterType;
typedef itk::MeshFileReader<MeshType> MeshReaderType;
typedef itk::MeanSquaresPointSetToImageMetric<PointSetType, DistanceImageType> MetricType;
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


struct Config {
    const static double imageMargin;
    static const double gradientTolerance;
    static const unsigned numberOfIterations;
};
const double Config::imageMargin = 50;
const double Config::gradientTolerance = 1e-5;
const unsigned Config::numberOfIterations = 100;



class Utils {
public:
    /**
     * read landmarks from the given file in slicer fcsv formant and return them as a list.
     *
     * The format is: label,x,y,z
     *
     * @param filename the filename
     * @returns A list of itk points
     */
    static std::vector<PointType > readLandmarks(const std::string& filename) {

        std::vector<PointType> ptList;

        std::fstream file ( filename.c_str() );
        if (!file) {
            std::cout << "could not read landmark file " << std::endl;
            throw std::runtime_error("could not read landmark file ");
        }
        std::string line;
        while (  std::getline ( file, line))
        {
            if (line.length() > 0 && line[0] == '#')
                continue;

            std::istringstream strstr(line);
            std::string token;
            std::getline(strstr, token, ','); // ignore the label
            std::getline(strstr, token, ','); // get the x coord
            double pt0 = atof(token.c_str());
            std::getline(strstr, token, ','); // get the y coord
            double pt1 = atof(token.c_str());
            std::getline(strstr, token, ','); // get the z coord
            double pt2 = atof(token.c_str());
            PointType pt;
            pt[0] = pt0; pt[1] = pt1; pt[2] = pt2;
            ptList.push_back(pt);
        }
        return ptList;
    }

};


// Returns a new model, that is restricted to go through the proints specified in targetLandmarks..
//
StatisticalModelType::Pointer
computePosteriorModel( const StatisticalModelType* statisticalModel,
                       const  std::vector<PointType >& modelLandmarks,
                      const  std::vector<PointType >& targetLandmarks,
                      double variance)
{

        // invert the transformand back transform the landmarks
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



//
// This class is used to track the progress of the optimization
// (its method Execute is called in each iteration of the optimization)
//
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


int main(int argc, char* argv[]) {


    if (argc < 4) {
        std::cout << "usage " << argv[0] << " modelname shape outputmesh [fixedLandmarks movingLandmarks lmVariance] " << std::endl;
        exit(-1);
    }

    char* modelName = argv[1];
    char* targetName = argv[2];
    char* outputMeshName = argv[3];

    char* fixedLandmarksName = 0;
    char* movingLandmarksName = 0;
    double lmVariance = 0;

    if (argc == 7) {
        fixedLandmarksName = argv[4];
        movingLandmarksName = argv[5];
        lmVariance = atof(argv[6]);
    }



    // load the Surface to which we will fit
    MeshReaderType::Pointer targetReader = MeshReaderType::New();
    targetReader->SetFileName(targetName);
    targetReader->Update();
    MeshType::Pointer mesh = targetReader->GetOutput();

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

    // Compute a bounding box around the reference shape
    typedef  itk::BoundingBox<int, 3, float, MeshType::PointsContainer> BoundingBoxType;
    BoundingBoxType::Pointer bb = BoundingBoxType::New();
    bb->SetPoints(representer->GetReference()->GetPoints());
    bb->ComputeBoundingBox();

    // Compute a binary image from the point set, which is as large as the bounding box plus a margin.
    PointsToImageFilterType::Pointer pointsToImageFilter = PointsToImageFilterType::New();
    pointsToImageFilter->SetInput( mesh );
    BinaryImageType::SpacingType spacing; spacing.Fill( 1.0 );
    BinaryImageType::PointType origin = bb->GetMinimum();
    BinaryImageType::SpacingType diff = bb->GetMaximum() - bb->GetMinimum();
    BinaryImageType::SizeType size;
    for (unsigned i =0; i < 3; i++) {
        origin[i] -= Config::imageMargin; // a five cm margin
        size[i] = diff[i] + Config::imageMargin;
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

    // Write out the fitting result
    itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
    writer->SetFileName(outputMeshName);
    writer->SetInput(transformMeshFilter->GetOutput());
    writer->Update();

}



