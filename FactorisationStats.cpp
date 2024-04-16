#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "FactorisationStats.h"


void FactorisationStats::OutputString(std::string s){
    #if extra_stats==1
    std::cout << s;
    if(filename != ""){
        std::ofstream outfile;
        outfile.open(filename, std::ios_base::app);
        outfile << s;
        outfile.close();
    }
    #endif
}

void FactorisationStats::OutputExtraStatsNumber(){
    #if extra_stats==1
    std::string stats_string = GetStatsForNumber();
    OutputString(stats_string);
    #endif
}

void FactorisationStats::OutputStatsAverage(){
    #if extra_stats==1
    std::string stats_string = GetStatsForAverage();
    OutputString(stats_string);
    #endif
}

void FactorisationStats::CreateFile(){
    #if extra_stats==1
    if(filename != ""){
        std::ofstream outfile;
        outfile.open(filename, std::ios_base::trunc);
        outfile << "START";
        outfile.close();
    }
    #endif
}

std::string FactorisationStats::GetStatsForNumber(){
    return "";
}

std::string FactorisationStats::GetStatsForAverage(){
    return "";
}

void StatsEcm::RunStart(double B1, double B2){
    #if extra_stats==1

    curves_per_thread = 0;
    total_time_s1 = 0;
    total_time_s2 = 0;
    
    this->B1 = B1;
    this->B2 = B2;
    #endif
}

void StatsEcm::RunEnd(){

}

void StatsEcm::AddCurveCount(int index){
    #if extra_stats==1
    curves_per_thread++;
    #endif
}

std::string StatsEcm::GetStatsForNumber(){
    #if extra_stats==1
    std::stringstream ss;

    int total_curves = 0;
    for(int i = 0; i < num_curves; i++){
        total_curves += curves_per_thread;
    }

    ss << "Total Curves: " << curves_per_thread;

    if(B2 != 0){
        double full_total_s1 = 0;
        double full_total_s2 = 0;

        full_total_s1 += total_time_s1;
        full_total_s2 += total_time_s2;


        if(full_total_s1 != 0 && full_total_s2 != 0){
        double K = (full_total_s1/B1)/(full_total_s2/(B2-B1));
        total_K += K;
        K_count++;

        ss << " Total time in S1: " << full_total_s1 << " Total time in S2: " << full_total_s2 << " K: " << K;
        }
    }

    nums_factored++;
    total_curves_per_size += total_curves;
    
    ss  << "\n";
    return ss.str();
    #else
    return "";
    #endif
}

std::string StatsEcm::GetStatsForAverage(){
    #if extra_stats==1
    std::stringstream ss;

    ss << "Avergaes Curve Count: " << (float)total_curves_per_size/nums_factored << " K: " << (float)total_K/K_count << "\n";

    total_K = 0;
    K_count = 0;
    nums_factored = 0;
    total_curves_per_size = 0;

    return ss.str();
    #else
    return "";
    #endif
}

void StatsEcm::StartTimingS1(int index){
    #if extra_stats==1
    s1_start = std::chrono::high_resolution_clock::now();
    #endif
}

void StatsEcm::EndTimingS1(int index){
    #if extra_stats==1
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end_time - s1_start;
    float timing = duration.count();
    curr_time_s1 = timing;
    #endif
}

void StatsEcm::StartTimingS2(int index){
    #if extra_stats==1
    s2_start = std::chrono::high_resolution_clock::now();
    #endif
}

void StatsEcm::EndTimingS2(int index){
    #if extra_stats==1
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = (end_time - s1_start);
    float timing = duration.count();

    total_time_s1 += curr_time_s1;
    total_time_s2 += timing;
    #endif
}
