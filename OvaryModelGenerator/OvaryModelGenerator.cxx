#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
//#include <vtkTransformFilter.h>
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
using namespace std;
int main(int argc, char* argv[])
{
  // verify number of parameters in command line
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " nslices nrot d h w outdir" << std::endl;
    std::cerr << "  nslices is an integer number of sections along the long-axis" << std::endl;
    std::cerr << "  nrot is an integer number of rotational sections (1, 2, or 4)" << std::endl;
    std::cerr << "  d is the anterior-posterior (shortest) ovary dimension in mm" << std::endl;
    std::cerr << "  h is the superior-inferior ovary dimension in mm" << std::endl; 
    std::cerr << "  w is the medial-lateral (longest) ovary dimension in mm" << std::endl; 
    std::cerr << "  outdir is the optional directory for the output files. Default is current directory" << std::endl;
    return EXIT_FAILURE;
    }

  // recover command line arguments
  int nslices  = atoi(argv[1]);
  int nrot = atoi(argv[2]);
  double d = atoi(argv[3]);
  double h = atoi(argv[4]);
  double w = atoi(argv[5]);

  const double scaleRatio = 0.65;

  // Experiment: shrink dims by 10
  double sd, sh, sw;
  sd = d * scaleRatio;
  sh = h * scaleRatio;
  sw = w * scaleRatio;

  const char *outdir = (argc < 7 ? "./output" : argv[6]);
  std::cout<<outdir<<endl;
  if (nslices != 1 && nslices != 3 && nslices != 12)
  {
    std:cerr<<"Error: nslices must be 1, 3, or 12"<<endl;
    return EXIT_FAILURE;
  }

  if (nrot != 1 && nrot != 2 && nrot != 4)
  {
    std::cerr<<"Error: Can only input 1, 2, or 4 in the nrot field"<<endl;
    return EXIT_FAILURE;
  }
  // if outdir does not include a '/', add it
  std::string outdirstr = outdir;
  outdirstr += ((outdir[strlen(outdir) - 1] == '/') ? "" : "/"); 

  std::string filename = "Ovary";
  std::string filename_p = "Pelvis";
  std::string fnimg = outdirstr + filename + ".nii";
  std::string fnmesh = outdirstr + filename + ".vtk";
  std::string fnxml = outdirstr + filename + ".vtp";
  std::string fnxmlpel = outdirstr + filename_p + ".vtp";
  std::string fnjson = outdirstr + filename;

  // reading in pelvis data

  //reads in data from the .glb file
  vtkNew<vtkGLTFReader> pelvisReader;
  pelvisReader->SetFileName("/Users/jileihao/data/data_hubmap/pelvis_model/VH_F_Pelvis.glb");
  pelvisReader->Update();

  //creates multiblock dataset to extract surface data from to convert into polydata
  vtkSmartPointer<vtkMultiBlockDataSet> pelvis_mb = pelvisReader->GetOutput(); 

  // specifying spheroid scale
  double Rsphere = 30;
  double sx = 0.5*sd/Rsphere;
  double sy = 0.5*sh/Rsphere;
  double sz = 0.5*sw/Rsphere;

  // create a spheroid with specified dimensions
  vtkNew<vtkSphereSource> sphereSource;
  sphereSource->SetThetaResolution(100);
  sphereSource->SetPhiResolution(100);
  sphereSource->SetRadius(Rsphere);
  vtkSmartPointer<vtkPolyData> pd = sphereSource->GetOutput();
  sphereSource->Update();

  //convert pelvis model to polydata
  //produces surface data from multiblock dataset
  vtkNew<vtkDataSetSurfaceFilter> pelvis_surfaces; 
  pelvis_surfaces->SetInputData(pelvis_mb);

  //compiles surface data into a single polydata object
  vtkNew<vtkCompositeDataGeometryFilter> pelvis_comp;
  pelvis_comp->SetInputConnection(pelvis_surfaces->GetOutputPort());
  pelvis_comp->Update();

  //rescaling parameters for ovary
  vtkNew<vtkTransform> rescale;
  rescale->Scale(sx, sy, sz);

  //rescale ovary polydata
  vtkNew<vtkTransformPolyDataFilter> rescaleFilter;
  rescaleFilter->SetInputConnection(sphereSource->GetOutputPort());
  rescaleFilter->SetTransform(rescale);
  rescaleFilter->Update();
  vtkSmartPointer<vtkPolyData> pd_trans = rescaleFilter->GetOutput();

  //rescale pelvis polydata
  vtkSmartPointer<vtkPolyData> pelvis_pd = pelvis_comp->GetOutput();

  double pel_bounds[6];
  pelvis_pd->GetBounds(pel_bounds);

  double x_width;

  double ovary_bounds[6];
  pd_trans->GetBounds(ovary_bounds);
  double z_ovaryDepth = ovary_bounds[5] - ovary_bounds[4];
  x_width = pel_bounds[1] - pel_bounds[0];

  //rescaling parameters for pelvis
  double pelScale = 290 / x_width; //average intercristal distance / original model width = scale factor to get a 29 cm (avg width) pelvis
  vtkNew<vtkTransform> pel_rescale;
  pel_rescale->Scale(pelScale, pelScale, pelScale);

  //rescale pelvis polydata
  vtkNew<vtkTransformPolyDataFilter> rescalePelvisFilter;
  rescalePelvisFilter->SetInputConnection(pelvis_comp->GetOutputPort());
  rescalePelvisFilter->SetTransform(pel_rescale);
  rescalePelvisFilter->Update();
  vtkSmartPointer<vtkPolyData> pelvisTransform_pd = rescalePelvisFilter->GetOutput();

  // create an image of the spheroid
  vtkNew<vtkImageData> whiteImage;  
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
    dim[i] = static_cast<int>(ceil((bounds[5] - bounds[4])/spacing[2])+1); 
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

        if (nrot == 1)
          {
          pix[0] = pix[0] + pix[0]*(slice_num-1);
          }

        if (nrot == 2)
          {        
          // determine which on side of the 45-deg plane the pixel is located 
          // and assign label
          double s45 = 0.7071*coords[0] + 0.7071*coords[1];
          if (s45 > 0)
            pix[0] = pix[0] + pix[0]*(slice_num-1); 
          else
            pix[0] = pix[0] + pix[0]*(slice_num-1+nslices);
          }

        if (nrot == 4)
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
  if (nrot == 0)
    imax = nslices;
  else
    imax = nrot * nslices;

  std::cout << "imax: " << imax << std::endl;
  std::vector<vtkSmartPointer<vtkPolyData>> vec;
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
    vec.push_back(labelMesh);
    }

  //Transforming pelvis model

  double baseMatrix[16] = {1, 0, 0, 0,
                           0, 1, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1}; 

  vtkTransform *pelTransform = vtkTransform::New();

  //Rotation
  pelTransform->RotateX(45);
  pelTransform->RotateY(-90);
  pelTransform->RotateZ(-10);
  
  //Translation
  pelTransform->Translate(-30, -35, 80);

  //Concatenate - allows for multiple transformation steps to be applied to a preset matrix
  pelTransform->Concatenate(baseMatrix);

  //Applying transformation to pelvis polydata
  vtkTransformPolyDataFilter *pelTransformFilter = vtkTransformPolyDataFilter::New();
  pelTransformFilter->SetTransform(pelTransform);
  pelTransformFilter->SetInputData(pelvisTransform_pd);
  pelTransformFilter->Update();

  //Retrieving output polydata of filter application
  vtkSmartPointer<vtkPolyData> pelFinal_pd = pelTransformFilter->GetOutput();

  // write vtp mesh
  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  ColorTable *ct = new ColorTable();
  vtkNew<vtkLookupTable> lut;
  lut->SetRange(1, imax);
  lut->SetSaturationRange(0.7, 1);
  lut->SetValueRange(0.7, 1);
  lut->Build();
  for(auto mesh : vec)
  {
      vtkNew<vtkActor> actor;
      vtkNew<vtkPolyDataMapper> mapper;
      mapper->SetInputData(mesh);

      auto fd = mesh->GetFieldData();
      auto label = fd->GetArray("Label")->GetTuple1(0);
      
      uint8_t r, g, b;
      float a;
      ct->GetLabelColor(label, r, g, b, a);

      double c[3];

      /*
      c[0] = r/255.0;
      c[1] = g/255.0;
      c[2] = b/255.0;
      */
      lut->GetColor(label, c);
      actor->GetProperty()->SetColor(c);
      actor->SetMapper(mapper);
      ren->AddActor(actor);
  }

  vtkNew<vtkActor> pvsActor;
  vtkNew<vtkPolyDataMapper> pvsMapper;
  pvsMapper->SetInputData(pelFinal_pd);
  pvsActor->SetMapper(pvsMapper);
  pvsActor->GetProperty()->SetOpacity(0.05);
  ren->AddActor(pvsActor);


  renWin->AddRenderer(ren);
  iren->SetRenderWindow(renWin);
  ren->SetBackground(0.32, 0.34, 0.43);
  ren->ResetCamera();
  iren->Initialize();
  
  
  // Display structure for troubleshooting
  /*
  vtkNew<vtkInteractorStyleSwitch> iSwitch;
  iSwitch->SetCurrentStyleToTrackballCamera();
  iren->SetInteractorStyle(iSwitch);
  ren->Render();
  iren->Start();
  */

  vtkNew<vtkJSONSceneExporter> exp;
  exp->SetRenderWindow(renWin);
  exp->SetActiveRenderer(ren);
  exp->SetFileName(fnjson.c_str());
  exp->Write();

  // Export a vtp scene for troubleshooting
  /*
  vtkNew<vtkSingleVTPExporter> vtpExp;
  vtpExp->SetRenderWindow(renWin);
  vtpExp->SetActiveRenderer(ren);
  std::string fnvtp = fnjson + ".vtp";
  vtpExp->SetFileName(fnvtp.c_str());
  vtpExp->Write();
  */

  return EXIT_SUCCESS;
}