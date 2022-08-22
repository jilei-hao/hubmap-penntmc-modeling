#ifndef __PelvisModelGenerator_h_
#define __PelvisModelGenerator_h_

#include <iostream>

#include "vtkPolyData.h"

class PelvisModelGenerator{

    public:
        float intercristalDistance;
        float pelvisScaleFactor;
        vtkPolyData* pelvisOutputPolydata;

    public:
        PelvisModelGenerator();
        ~PelvisModelGenerator();
        PelvisModelGenerator(const PelvisModelGenerator &other);

        void SetIntercristalDistancemm(float intercristalDistance); //sets intercristal distance of pelvis and scale factor
        void Generate();
        vtkPolyData* GetOutput();

};


#endif