#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

// Includes for pelvis model
  // First converted from .glb to a multiblock dataset
#include <vtkGLTFReader.h>
#include <vtkMultiBlockDataSet.h>

  // Then from multiblock data to surface data to polydata
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCompositeDataGeometryFilter.h>

#include "PelvisModelGenerator.h"

PelvisModelGenerator::PelvisModelGenerator(){

  this->intercristalDistance = 1;
  this->pelvisScaleFactor = 1;

};

PelvisModelGenerator::~PelvisModelGenerator(){};

PelvisModelGenerator::PelvisModelGenerator(const PelvisModelGenerator &other)
{

  this->intercristalDistance = other.intercristalDistance;
  this->pelvisScaleFactor = other.pelvisScaleFactor;
  this->pelvisOutputPolydata = other.pelvisOutputPolydata;

}

void PelvisModelGenerator::SetIntercristalDistancemm(float intercristalDistance){

  this->intercristalDistance = intercristalDistance;

};

void PelvisModelGenerator::Generate(){

    // reads in data from the .glb file
  vtkNew<vtkGLTFReader> pelvisReader;
  pelvisReader->SetFileName("./data/VH_F_Pelvis.glb");
  pelvisReader->Update();  

    // creates multiblock dataset to extract surface data from to convert into polydata
  vtkSmartPointer<vtkMultiBlockDataSet> pelvis_mb = pelvisReader->GetOutput(); 

  // convert pelvis model to polydata
    // produces surface data from multiblock dataset
  vtkNew<vtkDataSetSurfaceFilter> pelvis_surfaces; 
  pelvis_surfaces->SetInputData(pelvis_mb);

    // compiles surface data into a single polydata object
  vtkNew<vtkCompositeDataGeometryFilter> pelvis_comp;
  pelvis_comp->SetInputConnection(pelvis_surfaces->GetOutputPort());
  pelvis_comp->Update();

  // rescale pelvis polydata
    // retrieving pelvis polydata
  vtkSmartPointer<vtkPolyData> pelvis_pd = pelvis_comp->GetOutput();

    // retrieving bounds from that pelvis polydata to find modelWidth
  double pel_bounds[6];
  pelvis_pd->GetBounds(pel_bounds);
  double modelWidth = pel_bounds[1] - pel_bounds[0];
  
    // uses intercristalDistance from set function
    // setting scale factor
  this->pelvisScaleFactor = intercristalDistance / modelWidth; //average intercristal distance / original model width = scale factor to get a 29 cm (avg width) pelvis

  vtkNew<vtkTransform> pelvisRescale;
  pelvisRescale->Scale(pelvisScaleFactor, pelvisScaleFactor, pelvisScaleFactor);

    // apply rescale to filter
  vtkNew<vtkTransformPolyDataFilter> rescalePelvisFilter;
  rescalePelvisFilter->SetInputConnection(pelvis_comp->GetOutputPort());
  rescalePelvisFilter->SetTransform(pelvisRescale);
  rescalePelvisFilter->Update();
  vtkSmartPointer<vtkPolyData> pelvisTransform_pd = rescalePelvisFilter->GetOutput();


  // Transforming pelvis model

  double baseMatrix[16] = {1, 0, 0, 0,
                           0, 1, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1}; 

  vtkTransform *pelTransform = vtkTransform::New();

    // Rotation
  pelTransform->RotateX(45);
  pelTransform->RotateY(-90);
  pelTransform->RotateZ(-10);
  
    // Translation
  pelTransform->Translate(-30, -35, 80);

    // Concatenate - allows for multiple transformation steps to be applied to a preset matrix
  pelTransform->Concatenate(baseMatrix);

  // Applying transformation to pelvis polydata
  vtkTransformPolyDataFilter *pelTransformFilter = vtkTransformPolyDataFilter::New();
  pelTransformFilter->SetTransform(pelTransform);
  pelTransformFilter->SetInputData(pelvisTransform_pd);
  pelTransformFilter->Update();

  // Retrieving output polydata of filter application
  vtkSmartPointer<vtkPolyData> pelFinal_pd = pelTransformFilter->GetOutput();

  this->pelvisOutputPolydata = pelFinal_pd;

};

vtkPolyData* PelvisModelGenerator::GetOutput(){

  return this->pelvisOutputPolydata;

};