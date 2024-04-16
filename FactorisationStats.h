#pragma once
#include <string>
#include <chrono>

# define extra_stats 0 //collect extra stats off = 0 on = 1

class FactorisationStats{
    protected:
    std::string filename;
    
    public:
    FactorisationStats(std::string filename){
        this->filename = filename;
        CreateFile();
    }

    virtual std::string GetStatsForNumber();

    virtual std::string GetStatsForAverage();

    void CreateFile();

    void OutputExtraStatsNumber();

    void OutputStatsAverage();

    void OutputString(std::string s);
};

class StatsEcm : public FactorisationStats{
    public:

    #if extra_stats==1
    int curves_per_thread;
    int num_curves;
    int cumulative_curve_count;

    std::chrono::high_resolution_clock::time_point s1_start;
    std::chrono::high_resolution_clock::time_point s2_start;

    float curr_time_s1;
    float total_time_s1;
    float total_time_s2;

    double B1;
    double B2;
    double total_K;
    int total_curves_per_size;
    int nums_factored;
    int K_count;
    #endif

    StatsEcm(std::string filename, int num_curves) : FactorisationStats(filename){
        #if extra_stats==1
     //   curves_per_thread = new int[num_curves];
     //   s1_start = new std::chrono::high_resolution_clock::time_point[num_curves];
    //    s2_start = new std::chrono::high_resolution_clock::time_point[num_curves];
     //   curr_time_s1 = new float[num_curves];
    //    total_time_s1 = new float[num_curves];
     //   total_time_s2 = new float[num_curves];
        this->num_curves = num_curves;
        total_K = 0;
        nums_factored = 0;
        K_count = 0;
        total_curves_per_size = 0;
        #endif
    }

    ~StatsEcm(){
        #if extra_stats==1
      //  delete[] curves_per_thread;
      //  delete[] s1_start;
     //   delete[] s2_start;
     //   delete[] curr_time_s1;
     //   delete[] total_time_s1;
     //   delete[] total_time_s2;
        #endif
    }

    void AddCurveCount(int index);

    void RunStart(double B1, double B2);

    void RunEnd();

    std::string GetStatsForAverage();

    std::string GetStatsForNumber();

    void StartTimingS1(int index);

    void EndTimingS1(int index);

    void StartTimingS2(int index);

    void EndTimingS2(int index);
};

class StatsPollardsRho : public FactorisationStats{

};
