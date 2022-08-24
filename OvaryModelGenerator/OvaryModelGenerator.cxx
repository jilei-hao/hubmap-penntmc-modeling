#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkUnsignedShortArray.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include "TestingHelper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkNew.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkJSONSceneExporter.h"
#include "vtkSingleVTPExporter.h"
#include "vtkProperty.h"
#include "ColorTable.hxx"
#include "vtkLookupTable.h"
#include "OvaryModelGenerator.h"

// Includes for pelvis model
  //First converted from .glb to a multiblock dataset
#include <vtkGLTFReader.h>
#include <vtkMultiBlockDataSet.h>

  //Then from multiblock data to surface data to polydata
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCompositeDataGeometryFilter.h>
/**
 * This program generates an ellipsoid (closed surface, vtkPolyData) and converts it into volume
 * representation (vtkImageData) where the foreground voxels are 1 and the background voxels are
 * 0. The interior of the ellipsoid is then divided into nslice*nrot subsections: nrot long-axis slices, each
 * of which is rotationally subdivided into nrot subregions (1, 2, or 4). Each subregion is assigned a unique 
 * label. Internally vtkPolyDataToImageStencil is utilized. The resultant multi-label image is saved to 
 * disk in NiFTI file format (Ovary.nii) and vtkPolyData format (Ovary.vtk).
 */

OvaryModelGenerator::OvaryModelGenerator(int slice, int rot, double de, double he, double wi)
{
  this->SetLongAxisSlices(slice);
  this->SetRotationalSlices(rot);
  this->SetDimensions(de, he, wi);
}
OvaryModelGenerator::~OvaryModelGenerator(){};

void
OvaryModelGenerator::SetRotationalSlices(int rot)
{
  this->nrot = rot;
}

void
OvaryModelGenerator::SetLongAxisSlices(int slice)
{
  this->nslices = slice;
}
void
OvaryModelGenerator::SetDimensions(double de, double he, double wi)
{
  this->d = de * scaleRatio;
  this->h = he * scaleRatio;
  this->w = wi * scaleRatio;
}


void
OvaryModelGenerator::Generate()
{
  // compute dimensions
  double Rsphere = 30;
  double sx = 0.5*d/Rsphere;
  double sy = 0.5*h/Rsphere;
  double sz = 0.5*w/Rsphere;

  // create a spheroid with specified dimensions
  vtkSmartPointer<vtkSphereSource> sphereSource = 
    vtkSmartPointer<vtkSphereSource>::New();  
  sphereSource->SetThetaResolution(100);
  sphereSource->SetPhiResolution(100);
  sphereSource->SetRadius(Rsphere);
  vtkSmartPointer<vtkPolyData> pd = sphereSource->GetOutput();
  sphereSource->Update();

  vtkSmartPointer<vtkTransform> rescale = 
    vtkSmartPointer<vtkTransform>::New();
  rescale->Scale(sx, sy, sz);

  vtkSmartPointer<vtkTransformPolyDataFilter> rescaleFilter = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  rescaleFilter->SetInputConnection(sphereSource->GetOutputPort());
  rescaleFilter->SetTransform(rescale);
  rescaleFilter->Update();
  vtkSmartPointer<vtkPolyData> pd_trans = rescaleFilter->GetOutput();

  // create an image of the spheroid
  vtkSmartPointer<vtkImageData> whiteImage = 
    vtkSmartPointer<vtkImageData>::New();    
  double bounds[6];
  pd_trans->GetBounds(bounds);
  double spacing[3]; // desired volume spacing
  spacing[0] = 0.2;
  spacing[1] = 0.2;
  spacing[2] = 0.2;
  whiteImage->SetSpacing(spacing);

  // compute dimensions
  int dim[3];
  for (int i = 0; i < 3; i++)
    {
    //dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])+1);
    dim[i] = static_cast<int>(ceil((bounds[5] - bounds[4])/spacing[2])+1); 
    //std::cerr << "dim " << dim[i] << std::endl;
    //std::cerr << "bounds[] " << bounds[3] << " " << bounds[4] << " " << bounds[5] << endl;
    }
  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(-1, dim[0] - 1, -1, dim[1] - 1, -1, dim[2] - 1);

  double origin[3];
  origin[0] = bounds[4] + spacing[0] / 2;
  origin[1] = bounds[4] + spacing[1] / 2;
  origin[2] = bounds[4] + spacing[2] / 2;
  whiteImage->SetOrigin(origin);

  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);

  // fill the image with foreground voxels:
  unsigned char inval = 1;
  unsigned char outval = 0;
  vtkIdType count = whiteImage->GetNumberOfPoints();
  for (vtkIdType i = 0; i < count; ++i)
    {
    whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
    }

      // polygonal data --> image stencil:
      vtkNew<vtkPolyDataToImageStencil> pol2stenc;

  pol2stenc->SetInputData(pd_trans);
  pol2stenc->SetOutputOrigin(origin);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2stenc->Update();

  // cut the corresponding white image and set the background
  vtkNew<vtkImageStencil> imgstenc;

  imgstenc->SetInputData(whiteImage);
  imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue(outval);
  imgstenc->Update();

  // update pixel values along long axis of ellipsoid
  double step = (double)dim[2]/nslices;
  int slice_num = 1;

  // std::cout << "dim[2]=" << dim[2] << std::endl;
  // std::cout << "step=" << step << std::endl;

  vtkSmartPointer<vtkImageData> mlImage = imgstenc->GetOutput();

  for (int k = 0; k < dim[2]; k++)
    {
    if (k > step*slice_num)
      slice_num = slice_num + 1;
      
    for (int i = 0; i < dim[0]; i++)
      {
      for (int j = 0; j < dim[1]; j++)
        {
        // pixel value
        unsigned char* pix = 
          static_cast<unsigned char*>(mlImage->GetScalarPointer(i,j,k));

        // pixel physical coordinates 
        int ijk[3];
        ijk[0] = i;
        ijk[1] = j;
        ijk[2] = k;
        
        vtkIdType id = mlImage->ComputePointId(ijk);
        double* coords = 
          static_cast<double*>(mlImage->GetPoint(id));

        if (this->nrot == 1)
          {
          pix[0] = pix[0] + pix[0]*(slice_num-1);
          }

        if (this->nrot == 2)
          {        
          // determine which on side of the 45-deg plane the pixel is located 
          // and assign label
          double s45 = 0.7071*coords[0] + 0.7071*coords[1];
          if (s45 > 0)
            pix[0] = pix[0] + pix[0]*(slice_num-1); 
          else
            pix[0] = pix[0] + pix[0]*(slice_num-1+nslices);
          }

        if (this->nrot == 4)
          {
          double s45 = 0.7071*coords[0] + 0.7071*coords[1];
          double s135 = -0.7071*coords[0] + 0.7071*coords[1];
          if (s45 > 0)
            {
            if (s135 > 0)
              pix[0] = pix[0] + pix[0]*(slice_num-1+nslices);
            else
              pix[0] = pix[0] + pix[0]*(slice_num-1); 
            }
          else
            {
            if (s135 > 0)
              pix[0] = pix[0] + pix[0]*(slice_num-1+2*nslices);
            else
              pix[0] = pix[0] + pix[0]*(slice_num-1+3*nslices);
            }
          }

        }
      }
    }

/*
  // write the label image for troubleshooting
  vtkSmartPointer<vtkNIFTIImageWriter> writer = 
    vtkSmartPointer<vtkNIFTIImageWriter>::New();
  writer->SetFileName(fnimg.c_str());
  writer->SetInputData(mlImage);
  writer->Write();
*/


  // Run marching cubes on the image to convert it back to VTK polydata
  vtkPolyData *pipe_tail;

  // Extracting one label at a time and assigning label value
  float imax;
  if (this->nrot == 0)
    imax = this->nslices;
  else
    imax = this->nrot * this->nslices;

  std::cout << "imax: " << imax << std::endl;
  for (float i = 1; i <= imax; i += 1.0)
    {
    float lbl = floor(i);

    std::cout << "  -- Processing Label: " << lbl << std::endl;

    // Extract one label
    vtkNew<vtkDiscreteMarchingCubes> fltDMC;
    fltDMC->SetInputData(mlImage);
    fltDMC->ComputeGradientsOff();
    fltDMC->ComputeScalarsOff();
    fltDMC->SetNumberOfContours(1);
    fltDMC->ComputeNormalsOn();
    fltDMC->SetValue(0, lbl);
    fltDMC->Update();

    vtkSmartPointer<vtkPolyData> labelMesh = fltDMC->GetOutput();


    // Set scalar values for the label
    vtkNew<vtkUnsignedShortArray> scalar;
    scalar->SetNumberOfComponents(1);
    scalar->InsertNextTuple1(lbl);
    scalar->SetName("Label");
    labelMesh->GetFieldData()->AddArray(scalar);
    
    // Create the transform filter for ovary
    vtkTransformPolyDataFilter *fltTransform = vtkTransformPolyDataFilter::New();
    fltTransform->SetInputData(labelMesh);

    vtkSmartPointer<vtkTransform> vox2coords =
      vtkSmartPointer<vtkTransform>::New();
    vox2coords->Translate(0,0,0);

    // Experiment: Scale size back
    vtkNew<vtkTransform> scaleBack;
    double sbr = 1/scaleRatio; // scale back ratio
    scaleBack->Scale(sbr, sbr, sbr);
    vox2coords->Concatenate(scaleBack);

    // Update the VTK transform to match
    fltTransform->SetTransform(vox2coords);
    fltTransform->Update();
    labelMesh = fltTransform->GetOutput();

    // smoothing  
    
    vtkNew<vtkWindowedSincPolyDataFilter> smoothFilter;
    smoothFilter->SetInputData(labelMesh);
    smoothFilter->SetNumberOfIterations(50);
    smoothFilter->SetPassBand(0.05);
    smoothFilter->SetNonManifoldSmoothing(false);
    smoothFilter->Update();
    labelMesh = smoothFilter->GetOutput();
    
    // Get final output
    m_Output.push_back(labelMesh);
    }
}

std::vector<vtkSmartPointer<vtkPolyData>> &
OvaryModelGenerator
::GetOutput()
{
  return m_Output;
}

// void
// OvaryModelGenerator::GetOutput(std::vector<vtkSmartPointer<vtkPolyData>> vec, int imax)
// {
//   // vec = Generate();

//   
  
  
//   // Display structure for troubleshooting
//   /*
//   vtkNew<vtkInteractorStyleSwitch> iSwitch;
//   iSwitch->SetCurrentStyleToTrackballCamera();
//   iren->SetInteractorStyle(iSwitch);
//   ren->Render();
//   iren->Start();
//   */
// }