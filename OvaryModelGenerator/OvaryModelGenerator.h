#ifndef __OvaryModelGenerator_h_
#define __OvaryModelGenerator_h_

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <vector>

class OvaryModelGenerator
{
    public:
        OvaryModelGenerator(int slice, int rot, double de, double he, double wi);
        ~OvaryModelGenerator();
        OvaryModelGenerator(const OvaryModelGenerator &other) = delete;
        OvaryModelGenerator & operator= (const OvaryModelGenerator & other) = delete;
        // static int inputCheck();
        void SetRotationalSlices(int rot);
        void SetOutDir(std::string out);
        void SetLongAxisSlices(int slice);
        void SetDimensions(double de, double he, double wi);
        void Generate();
        std::vector<vtkSmartPointer<vtkPolyData>> &GetOutput();

    private:
        int nslices, nrot;
        double d, h, w;
        const double scaleRatio = 0.65;
        std::vector<vtkSmartPointer<vtkPolyData>> m_Output;

};

#endif