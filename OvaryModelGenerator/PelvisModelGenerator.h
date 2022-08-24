#ifndef __PelvisModelGenerator_h_
#define __PelvisModelGenerator_h_

#include <iostream>

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

class PelvisModelGenerator{

    public:
        PelvisModelGenerator();
        ~PelvisModelGenerator();
        PelvisModelGenerator(const PelvisModelGenerator &other) = delete;
        PelvisModelGenerator & operator=(const PelvisModelGenerator &other) = delete;

        void SetIntercristalDistancemm(float intercristalDistance); //sets intercristal distance of pelvis and scale factor
        void Generate();
        vtkPolyData* GetOutput();

    private:
        float m_IntercristalDistance;
        float m_PelvisScaleFactor;
        vtkSmartPointer<vtkPolyData> m_PelvisOutputPolydata;
};


#endif