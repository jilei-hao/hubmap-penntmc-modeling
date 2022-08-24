#include <iostream>
#include <sstream>

#include "OvaryModelGenerator.h"
#include "PelvisModelGenerator.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkNew.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkJSONSceneExporter.h"
#include "vtkSingleVTPExporter.h"
#include "vtkFieldData.h"


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
  const char *outdir = (argc < 7 ? "./output" : argv[6]);
  std::string outdirstr = outdir;
  outdirstr += ((outdir[strlen(outdir) - 1] == '/') ? "" : "/"); 
  // recover command line arguments
  int nslices  = atoi(argv[1]);
  int nrot = atoi(argv[2]);
  double d = atoi(argv[3]);
  double h = atoi(argv[4]);
  double w = atoi(argv[5]);
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

    // Generate pelvis model
  const float pelvisIntercristalDistancemm = 290;
    
  PelvisModelGenerator *pvsGen = new PelvisModelGenerator();
  pvsGen->SetIntercristalDistancemm(pelvisIntercristalDistancemm);
  pvsGen->Generate();
  vtkSmartPointer<vtkPolyData> pelFinal_pd = pvsGen->GetOutput();
    
  delete pvsGen;

  OvaryModelGenerator myOvary(nslices, nrot, d, h, w);
  myOvary.Generate();
  auto ovarySections = myOvary.GetOutput();

  // Model Writing Logic
  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  for(auto mesh : ovarySections)
  {
    vtkNew<vtkActor> actor;
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(mesh);
    actor->SetMapper(mapper);
    ren->AddActor(actor);
  }
  renWin->AddRenderer(ren);
  iren->SetRenderWindow(renWin);
  ren->SetBackground(0.32, 0.34, 0.43);
  ren->ResetCamera();
  iren->Initialize();

  std::ostringstream fnOut;
  fnOut << outdirstr << "Ovary";

  vtkNew<vtkJSONSceneExporter> exporter;
  exporter->SetInput(renWin);
  exporter->SetActiveRenderer(ren);
  exporter->SetFileName(fnOut.str().c_str());
  exporter->Write();

  return EXIT_SUCCESS;
}
