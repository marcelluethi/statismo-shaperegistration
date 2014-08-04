/*
 *
 * Author: Marcel Luethi, (marcel.luethi@unibas.ch)
 *
 * Copyright (c) 2011 University of Basel
 * All rights reserved.
 *
 */

// If you have a c++11 compatible compiler, you can uncomment the following line.
// This will parallelize the computation and thus greatly improve performance of the model building.
// (This also requires that ITK was compiled with C++11 support)
//#define HAS_CXX11_ASYNC 1

#include "itkStandardMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkLowRankGPModelBuilder.h"
#include "itkMeshFileReader.h"
#include "itkMesh.h"




/**
 * The kernel defines our prior assumptions about the deformations in the registration.
 * Here we use a simple Gaussian Kernel, which leads to smooth deformations. We assume that the
 * components of the deformation field are uncorrelated.
 *
 * The variance of the kernel defines how smooth the deformations are and the scale the size (in mm)
 * of an average deformation.
*/
template <class TPoint>
class GaussianKernel: public statismo::MatrixValuedKernel<TPoint>{
public:
    typedef typename  TPoint::CoordRepType CoordRepType;
    typedef vnl_vector<CoordRepType> VectorType;

    GaussianKernel(double variance, double scale) : statismo::MatrixValuedKernel<TPoint>(3), m_variance(variance), m_scale(scale) {    }

    inline statismo::MatrixType operator()(const TPoint& x, const TPoint& y) const {
        VectorType xv = x.GetVnlVector();
        VectorType yv = y.GetVnlVector();

        VectorType r = yv - xv;
        return statismo::MatrixType::Identity(3,3) * m_scale * exp(-dot_product(r, r) / m_variance);
    }

    std::string GetKernelInfo() const {
        std::ostringstream os;
        os << "GaussianKernel(" << m_variance << ")" << "* " << m_scale;
        return os.str();
    }

private:

    double m_variance;
    double m_scale;
};




int main(int argc, char *argv[]) {

    typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
    typedef RepresenterType::DatasetType MeshType;
    typedef typename RepresenterType::PointType PointType;

    typedef itk::LowRankGPModelBuilder<MeshType> ModelBuilderType;
    typedef itk::StatisticalModel<MeshType> StatisticalModelType;
    typedef itk::MeshFileReader<MeshType>  MeshFileReaderType;

    const unsigned numNystromPoints = 500;

    if (argc != 6) {
        std::cout << "usage\t" << argv[0] << " referenceFilename kernelWidth kernelScale numberOfBasisFunctions outputFilename" << std::endl;
        exit(-1);
    }

    std::string referenceFn = argv[1];
    double kernelWidth = std::atof(argv[2]);
    double kernelScale = std::atof(argv[3]);
    double numBasisFunctions = std::atoi(argv[4]);
    std::string outputFn = argv[5];

    typename MeshFileReaderType::Pointer referenceReader = MeshFileReaderType::New();
    referenceReader->SetFileName(referenceFn);
    referenceReader->Update();

    typename RepresenterType::Pointer representer = RepresenterType::New();
    representer->SetReference(referenceReader->GetOutput());

    const GaussianKernel<PointType> gk = GaussianKernel<PointType>(kernelWidth * kernelWidth, kernelScale);

    typename ModelBuilderType::Pointer gpModelBuilder = ModelBuilderType::New();
    gpModelBuilder->SetRepresenter(representer);
    typename StatisticalModelType::Pointer model = gpModelBuilder->BuildNewModel(representer->GetReference(), gk, numBasisFunctions, numNystromPoints);
    model->Save(outputFn.c_str());
}

