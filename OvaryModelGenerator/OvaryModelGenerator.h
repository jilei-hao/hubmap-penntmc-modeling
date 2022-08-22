#ifndef __OvaryModelGenerator_h_
#define __OvaryModelGenerator_h_

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <vector>

class OvaryModelGenerator
{
    public:
        int nslices, nrot;
        double d, h, w;
        std::string outdirstr;

        OvaryModelGenerator(int slice, int rot, double de, double he, double wi, std::string outdir);
        ~OvaryModelGenerator();
        OvaryModelGenerator(const OvaryModelGenerator &other);
        // static int inputCheck();
        void SetRotationalSlices(int rot);
        void SetOutDir(std::string out);
        void SetLongAxisSlices(int slice);
        void SetDimensions(double de, double he, double wi);
        void Generate();
        void GetOutput(std::vector<vtkSmartPointer<vtkPolyData>> vec, int imax);







};

#endif