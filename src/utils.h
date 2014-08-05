/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
  *
 * Copyright (c) 2014 University of Basel
 * All rights reserved.
 */


#ifndef UTILS_H
#define UTILS_H

#include "itkMeshFileWriter.h"
#include "itkMeshFileReader.h"
#include "itkIdentityTransform.h"
#include "itkTransformMeshFilter.h"

typedef itk::Point<double, 3> PointType;

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


    template<typename MeshType>
    static void writeMesh(MeshType* mesh, const std::string& outputMeshName) {
      typename itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
       writer->SetFileName(outputMeshName);
       writer->SetInput(mesh);
       writer->Update();

    }

    template<typename MeshType>
    static typename MeshType::Pointer readMesh(const std::string& filename) {
        typename itk::MeshFileReader<MeshType>::Pointer reader = itk::MeshFileReader<MeshType>::New();
        reader->SetFileName(filename);
        reader->Update();
        typename MeshType::Pointer mesh = reader->GetOutput();
        return mesh;
    }


    // clones a mesh
    template<typename MeshType>
    static typename MeshType::Pointer cloneMesh(MeshType* mesh) {

        // as there is no clone method in itk, we use the transform filter with an identity transform,
        // which effectively clones the mesh.
        typedef itk::IdentityTransform<float, MeshType::PointDimension> IdentityTransformType;
        typedef itk::TransformMeshFilter<MeshType, MeshType, IdentityTransformType> TransformFilterType;

        typename TransformFilterType::Pointer transformMeshFilter = TransformFilterType::New();
        typename IdentityTransformType::Pointer id = IdentityTransformType::New();
        transformMeshFilter->SetInput(mesh);
        transformMeshFilter->SetTransform(id);
        transformMeshFilter->Update();
        return transformMeshFilter->GetOutput();

    }


};


#endif // UTILS_H
