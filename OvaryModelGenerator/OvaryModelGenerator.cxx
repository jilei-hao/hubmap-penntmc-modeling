#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>

/**
 * This program generates an ellipsoid (closed surface, vtkPolyData) and converts it into volume
 * representation (vtkImageData) where the foreground voxels are 1 and the background voxels are
 * 0. The interior of the ellipsoid is then divded into 48 subsections: 12 long-axis slices, each
 * of which is rotationally subdivided into four subregions.  
 * Internally vtkPolyDataToImageStencil is utilized. The resultant multi-label image is saved to 
 * disk in NiFTI file format (Ovary.nii).
 */
int main(int, char *[])
{  
  vtkSmartPointer<vtkSphereSource> sphereSource = 
    vtkSmartPointer<vtkSphereSource>::New();  
  sphereSource->SetRadius(20);
  sphereSource->SetPhiResolution(50);
  sphereSource->SetThetaResolution(50);
  vtkSmartPointer<vtkPolyData> pd = sphereSource->GetOutput();
  sphereSource->Update();

  vtkSmartPointer<vtkTransform> rescale = 
    vtkSmartPointer<vtkTransform>::New();
  rescale->Scale(0.5, 0.8, 1.2);

  vtkSmartPointer<vtkTransformPolyDataFilter> rescaleFilter = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  rescaleFilter->SetInputConnection(sphereSource->GetOutputPort());
  rescaleFilter->SetTransform(rescale);
  rescaleFilter->Update();
  vtkSmartPointer<vtkPolyData> pd_trans = rescaleFilter->GetOutput();

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
   }
  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(-1, dim[0] - 1, -1, dim[1] - 1, -1, dim[2] - 1);

  double origin[3];
  origin[0] = bounds[4] + spacing[0] / 2;
  origin[1] = bounds[4] + spacing[1] / 2;
  origin[2] = bounds[4] + spacing[2] / 2;
  whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
  whiteImage->SetScalarTypeToUnsignedChar();
  whiteImage->AllocateScalars();
#else
  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
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
#if VTK_MAJOR_VERSION <= 5
  pol2stenc->SetInput(pd_trans);
#else
  pol2stenc->SetInputData(pd_trans);
#endif
  pol2stenc->SetOutputOrigin(origin);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2stenc->Update();

  // cut the corresponding white image and set the background:
  vtkSmartPointer<vtkImageStencil> imgstenc = 
    vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
  imgstenc->SetInput(whiteImage);
  imgstenc->SetStencil(pol2stenc->GetOutput());
#else
  imgstenc->SetInputData(whiteImage);
  imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
  imgstenc->ReverseStencilOff();
  imgstenc->SetBackgroundValue(outval);
  imgstenc->Update();

  // update pixel values along long axis of ellipsoid
  int step = dim[2]/12;
  int slice_num = 1;
  vtkSmartPointer<vtkImageData> mlImage = imgstenc->GetOutput();
  for (int k = 0; k < dim[2]; k++)
    {

      if (k >= step*slice_num)
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
        
        // determine which on side of the 45-deg plane the pixel is located 
        // and assign label
        double s45 = 0.7071*coords[0] + 0.7071*coords[1];
        double s135 = -0.7071*coords[0] + 0.7071*coords[1];
        if (s45 > 0){
            if (s135 > 0)
                pix[0] = pix[0] + pix[0]*(slice_num-1+12);
            else
                pix[0] = pix[0] + pix[0]*(slice_num-1); //
            }
        else{
            if (s135 > 0)
                pix[0] = pix[0] + pix[0]*(slice_num-1+2*12);
            else
                pix[0] = pix[0] + pix[0]*(slice_num-1+3*12);
            }
            }
     	}
    }

  // write the label image
  vtkSmartPointer<vtkNIFTIImageWriter> writer = 
    vtkSmartPointer<vtkNIFTIImageWriter>::New();
  writer->SetFileName("Ovary.nii");
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(mlImage);
#else
  writer->SetInputData(mlImage);
#endif
  writer->Write();  
  
  return EXIT_SUCCESS;
}
