#ifndef __TestingHelper_h_
#define __TestingHelper_h_

#include <iostream>
#include <vector>

#include "vtkPolyData.h"

class TestingHelper
{
public:
    TestingHelper();
    ~TestingHelper();

    /** 
     * Create a render window that renders the meshes in the list
     * Example:
     *   std::vector<vtkPolyData*> meshes;
     *   meshes.push_back(YourMeshData);
     *   TestingHelper::RenderMeshes(meshes);
     * 
    */
    static void RenderMeshes(std::vector<vtkPolyData*> &meshes);
};


#endif