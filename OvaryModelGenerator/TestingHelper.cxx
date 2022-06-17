#include "TestingHelper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkNew.h"
#include "vtkInteractorStyleSwitch.h"

TestingHelper
::TestingHelper()
{
    std::cout << "[TestingHelper] constructor" << std::endl;
}

TestingHelper
::~TestingHelper()
{
    std::cout << "[TestingHelper] destructor" << std::endl;
}

void
TestingHelper
::RenderMeshes(std::vector<vtkPolyData*> &meshes)
{
    std::cout << "[TestingHelper] RenderMeshes. Number of meshes=" << meshes.size() << std::endl;

    vtkNew<vtkRenderer> ren;
    vtkNew<vtkRenderWindow> renWin;
    vtkNew<vtkRenderWindowInteractor> iren;

    
    for(auto mesh : meshes)
    {
        vtkNew<vtkActor> actor;
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(mesh);
        actor->SetMapper(mapper);
        ren->AddActor(actor);
    }
    
    renWin->AddRenderer(ren);

    iren->SetRenderWindow(renWin);
    vtkNew<vtkInteractorStyleSwitch> iSwitch;
    iSwitch->SetCurrentStyleToTrackballCamera();
    iren->SetInteractorStyle(iSwitch);

    renWin->SetSize(800, 600);
    ren->SetBackground(0.4, 0.5, 0.6);
    ren->ResetCamera();
    renWin->Render();
    
    iren->Initialize();
    iren->Start();
}