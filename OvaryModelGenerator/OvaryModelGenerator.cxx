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

/**
 * This program generates an ellipsoid (closed surface, vtkPolyData) and converts it into volume
 * representation (vtkImageData) where the foreground voxels are 1 and the background voxels are
 * 0. The interior of the ellipsoid is then divided into nslice*nrot subsections: nrot long-axis slices, each
 * of which is rotationally subdivided into nrot subregions (1, 2, or 4). Each subregion is assigned a unique 
 * label. Internally vtkPolyDataToImageStencil is utilized. The resultant multi-label image is saved to 
 * disk in NiFTI file format (Ovary.nii) and vtkPolyData format (Ovary.vtk).
 */

int main(int argc, char* argv[])
{
  // verify number of parameters in command line
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " nslices nrot d h w fnimg fnmesh" << std::endl;
    std::cerr << "  nslices is an integer number of sections along the long-axis" << std::endl;
    std::cerr << "  nrot is an integer number of rotational sections (1, 2, or 4)" << std::endl;
    std::cerr << "  d is the anterior-posterior (shortest) ovary dimension in mm" << std::endl;
    std::cerr << "  h is the superior-inferior ovary dimension in mm" << std::endl; 
    std::cerr << "  w is the medial-lateral (longest) ovary dimension in mm" << std::endl; 
    std::cerr << "  outdir is the directory for the output files. Default is current directory" << std::endl;
    return EXIT_FAILURE;
    }

  // recover command line arguments
  int nslices  = atoi(argv[1]);
  int nrot = atoi(argv[2]);
  double d = atoi(argv[3]);
  double h = atoi(argv[4]);
  double w = atoi(argv[5]);
  const char *outdir = (argc < 7 ? "./output" : argv[6]);

  // if outdir does not include a '/', add it
  std::string outdirstr = outdir;
  outdirstr += ((outdir[strlen(outdir) - 1] == '/') ? "" : "/"); 

  std::string filename = "Ovary";
  std::string fnimg = outdirstr + filename + ".nii";
  std::string fnmesh = outdirstr + filename + ".vtk";
  std::string fnxml = outdirstr + filename + ".vtp";
  std::string fnjson = outdirstr + filename;

  //std::cout << "fnimg=" << fnimg << std::endl;
  //std::cout << "fnmesh=" << fnmesh << std::endl;

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
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = 
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();

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



  // write the label image for troubleshooting
  vtkSmartPointer<vtkNIFTIImageWriter> writer = 
    vtkSmartPointer<vtkNIFTIImageWriter>::New();
  writer->SetFileName(fnimg.c_str());
  writer->SetInputData(mlImage);
  writer->Write();


  // Run marching cubes on the image to convert it back to VTK polydata
  vtkPolyData *pipe_tail;

  // Append filter for assembling labels
  vtkAppendPolyData *fltAppend = vtkAppendPolyData::New();

  // Extracting one label at a time and assigning label value
  float imax;
  if (nrot == 0)
    imax = nslices;
  else
    imax = nrot * nslices;

  std::cout << "imax: " << imax << std::endl;

  for (float i = 1; i <= imax; i += 1.0)
    {
    
    float lbl = floor(i);

    std::cout << "  -- Processing Label: " << lbl << std::endl;

    // Extract one label
    vtkDiscreteMarchingCubes *fltDMC = vtkDiscreteMarchingCubes::New();
    fltDMC->SetInputData(mlImage);
    fltDMC->ComputeGradientsOff();
    fltDMC->ComputeScalarsOff();
    fltDMC->SetNumberOfContours(1);
    fltDMC->ComputeNormalsOn();
    fltDMC->SetValue(0, lbl);
    fltDMC->Update();

    vtkPolyData *labelMesh = fltDMC->GetOutput();

    // Set scalar values for the label
    vtkUnsignedShortArray *scalar = vtkUnsignedShortArray::New();
    scalar->SetNumberOfComponents(1);
    for (vtkIdType i = 0; i < labelMesh->GetNumberOfPoints(); ++i)
      {
      scalar->InsertNextTuple1(lbl);
      }

    scalar->SetName("Label");
    labelMesh->GetPointData()->SetScalars(scalar);
    fltAppend->AddInputData(labelMesh);
    }

  fltAppend->Update();
  pipe_tail = fltAppend->GetOutput();

  // Create the transform filter
  vtkTransformPolyDataFilter *fltTransform = vtkTransformPolyDataFilter::New();
  fltTransform->SetInputData(pipe_tail);
 
  // Compute the transform from VTK coordinates to NIFTI/RAS coordinates
  // typedef vnl_matrix_fixed<double, 4, 4> Mat44;
  // Mat44 vtk2out = ConstructVTKtoNiftiTransform(
  // mlImage->GetDirection().GetVnlMatrix(),
  // mlImage->GetOrigin().GetVnlVector(),
  // mlImage->GetSpacing().GetVnlVector());

  vtkSmartPointer<vtkTransform> vox2coords =
    vtkSmartPointer<vtkTransform>::New();
  //vox2coords->Scale(spacing[0],spacing[1],spacing[2]);
  vox2coords->Translate(0.2,0.2,0.2);

  // Update the VTK transform to match
  // vtkTransform *transform = vtkTransform::New();
  // transform->SetMatrix(vtk2out.data_block());
  fltTransform->SetTransform(vox2coords);
  fltTransform->Update();

  // Get final output
  vtkPolyData *mesh = fltTransform->GetOutput();

  // Flip normals if determinant of SFORM is negative
  //if(transform->GetMatrix()->Determinant() < 0)
  //  {
  //  vtkPointData *pd = mesh->GetPointData();
  //  vtkDataArray *nrm = pd->GetNormals();
  //  for(size_t i = 0; i < (size_t)nrm->GetNumberOfTuples(); i++)
  //    for(size_t j = 0; j < (size_t)nrm->GetNumberOfComponents(); j++)
  //      nrm->SetComponent(i,j,-nrm->GetComponent(i,j));
  //  nrm->Modified();
  //  }

  // write the label mesh
  /*
  vtkSmartPointer<vtkPolyDataWriter> meshWriter = 
    vtkSmartPointer<vtkPolyDataWriter>::New();
  
  meshWriter->SetFileName(fnmesh.c_str());
  meshWriter->SetInputData(mesh);
  meshWriter->Write();
  */

  // smoothing
  vtkNew<vtkWindowedSincPolyDataFilter> smoothFilter;
  smoothFilter->SetInputData(mesh);
  smoothFilter->SetNumberOfIterations(20);
  smoothFilter->SetPassBand(0.001);
  smoothFilter->SetNonManifoldSmoothing(false);
  smoothFilter->Update();
  
  // write vtp mesh
  vtkNew<vtkXMLPolyDataWriter> xmlWriter;
  xmlWriter->SetFileName(fnxml.c_str());
  xmlWriter->SetInputData(smoothFilter->GetOutput());
  xmlWriter->Write();

  return EXIT_SUCCESS;
}
