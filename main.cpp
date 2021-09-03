#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <bitset>
#include "time.h"
#include <cmath>

const std::vector<std::string> FILES{"Distances.dat"};
const int MAX_PROBLEM_SIZE = 100;
typedef std::bitset<MAX_PROBLEM_SIZE> Solution;


#define GRASP_TIME_LIMIT 10;
#define GRASP_ITER_LIMIT 10000000000;
#define ALPHA 0.80;
#define JUMP_LEVEL 2;

enum NEIGHBOURHOOD_OPERATOR {
    FLIP,
    DOUBLE_FLIP,
    TRIPLE_FLIP,
    QUAD_FLIP,
    DOUBLE_SWAP,
    DOUBLE_SWAP_SWAP,
    DOUBLE_FLIP_SWAP,
    SWAP
};

struct Candidate    {
    Candidate(): i(-1), contribution(0) {}
    Candidate(int j, int w): i(j), contribution(w) {}
    int i;
    int contribution;
};



void printSolution(int n, Solution x)  {

    printf("The Maximum Diversity Solution is to choose candidates ");
    for (int i = x._Find_first(); i < n; i = x._Find_next(i)) {
        printf("%d, ", i);
    }
    printf("\n");

}


//read function from lectures
bool readInstance(std::string &filename, std::vector<std::vector<int> >& a, int &n, int &m) {
    bool success = true;
    a.clear();
    // open the file:
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        // failed to open file, so make sure to indicate this when we return:
        success = false;
    } else {
        // file is open and we are ready to read:
        int temp;
        file >> n >> m;

        a.resize(n);

        for (int i = 0; i < n; ++i) {

            for (int j = 0; j < n; ++j) {
                file >> temp;
                a[i].push_back(temp);
            }
        }
    }
    return success;
}
//
//bool checkSolutionFeasibility(const int&n, const int&m, const Solution& x, std::vector<std::vector<int>>& a,
//                              std::vector<int>& b)  {
//    std::vector<int> row_activity(m);
//    computeSolutionRowActivities(n, m, x, a, b, row_activity);
//    for (int i = 0; i < m; ++i)
//        if (row_activity[i] > b[i])
//            return false;
//    return true;
//}


bool compCandidates(Candidate* a, Candidate* b) {
    return a->contribution < b->contribution;
}


void createRestrictedCandidateList(const int n, const std::vector<Candidate*>& candidateList, Solution& x, double alpha,
                                   std::vector<Candidate*>& restrictedCandidateList) {

    auto maxElement = std::max_element(candidateList.begin(), candidateList.end(), compCandidates);
    int contMax = ((*maxElement)->contribution);
    int k = 0;
    for (int i = 0; i < n; ++i) {
        if (candidateList.at(i)->contribution >= alpha * contMax) {
            restrictedCandidateList.push_back(candidateList.at(i));
            k++;
        }
    }
    int rclsize = restrictedCandidateList.size();
    std::sort(restrictedCandidateList.begin(), restrictedCandidateList.end(), compCandidates);
}

void freeCandidateMemory(std::vector<Candidate*> candidateList, int n)  {
    for (int i = 0; i < n; ++i)
        delete candidateList.at(i);
}

double constructionHeuristic(const int n, const int m, Solution& x, int& objVal, double alpha, std::vector<std::vector<int>>& a) {

    clock_t start = clock(); //start time

    x.reset();
    srand(time(NULL));
    int chosen_cand = rand() % n;            // choose first one randomly
    x.set(chosen_cand);

    std::vector<Candidate*> candidateList;
    objVal = 0;

    for (int i = 0; i < n; ++i)
        candidateList.push_back(new Candidate(i, a[i][chosen_cand]));

    //start algorithm
    for (int new_candidate_index = 1; new_candidate_index < m; ++new_candidate_index) {

        std::vector<Candidate*> restrictedCandidateList;
        createRestrictedCandidateList(n, candidateList, x, alpha, restrictedCandidateList);
        auto RCL_size = restrictedCandidateList.size();

        Candidate* new_candidate = restrictedCandidateList[rand() % RCL_size];
        x.set(new_candidate->i);
        objVal += new_candidate->contribution;
        for (int i = 0; i < n; ++i)
            candidateList.at(i)->contribution += a[i][new_candidate->i];
        new_candidate->contribution = 0;
    }


    freeCandidateMemory(candidateList, n);

    return (clock() - start) / double(CLOCKS_PER_SEC / 1000 ); //time in milliseconds spent on processing
    // (now time - start time)
}

int computeSolutionObjective(const int &n, std::vector<std::vector<int>>& a, const Solution& solution) {
    int objectiveValue = 0;
    for (auto i = solution._Find_first(); i < n; i = solution._Find_next(i))  {
        for (auto j = solution._Find_first(); j < n; j = solution._Find_next(j))  {
            objectiveValue += a[i][j];
        }
    }
    return objectiveValue/2;
}




bool yieldNeighbour(const int &n, const int &m,  Solution& x, int &objVal, std::vector<std::vector<int>>& a,
                    int jump_level = 2) {
#ifdef DEBUG_INFO
    freopen("log.out","a+", stdout);
#endif



    int global_improvement = 0;
    int best_out_gl;
    int best_in_gl;
    int best_out_fl;
    int best_in_fl;

    for (auto i = x._Find_first(); i < n; i = x._Find_next(i)) {

        int objectiveDecrease = 0;

        for (auto j = x._Find_first(); j < n; j = x._Find_next(j))
            objectiveDecrease += a[i][j];

        int first_level_improvement = 0;

        for (int j = 0; j < n; ++j) {

            if (x.test(j))
                continue;
            int objImprovement;
            int objectiveIncrease = 0;
            for (auto k = x._Find_first();
                 k < n; k = x._Find_next(k)) {    // compute what the entering value may give to the solution
                if (k != i)
                    objectiveIncrease += a[k][j];
            }

            if ((objImprovement = objectiveIncrease - objectiveDecrease) > 0) {
                if (objImprovement > global_improvement) {
                    global_improvement = objImprovement;
                    best_out_gl = i;
                    best_in_gl = j;
                }
                if (objImprovement > first_level_improvement) {
                    first_level_improvement = objImprovement;
                    best_out_fl = i;
                    best_in_fl = j;
                }



                if (jump_level == 2) {
                    x.set(i, 0);
                    x.set(j, 1);
                    objVal += objImprovement;
                    return true;
                }
            }
        }

        if (first_level_improvement > 0 && jump_level == 1) {
            x.set(best_in_fl, 1);
            x.set(best_out_fl, 0);
            objVal += first_level_improvement;
            return true;
        }

    }

    if (global_improvement > 0 && jump_level == 0) {
        x.set(best_in_gl, 1);
        x.set(best_out_gl, 0);
        objVal += global_improvement;
        return true;
    }


    return false;

}

double localSearch(const int &n, const int &m,  Solution& x, int &objVal,
                 std::vector<std::vector<int>>& a, bool& improved, double timeToFinish,
                 int jump_level = 2) {

    clock_t start = clock();
    improved = false;
    while (yieldNeighbour(n, m, x, objVal, a, jump_level))
    {
        improved = true;
        if (clock() > timeToFinish)
            break;
    }

    auto time_taken = (clock() - start) / double(CLOCKS_PER_SEC / 1000); // time taken by the heuristic in milliseconds

    return time_taken;
}


long long GRASP(const int &n, const int &m, Solution& x,
           std::vector<std::vector<int>> &a, double alpha, int  jumplevel = 2,
           int grasp_time_limit = 1000, int grasp_iter_limit = 100000) {


    std::vector<int> restrictedCandList(n);

    auto timeToFinish = double(CLOCKS_PER_SEC)*grasp_time_limit + clock();

    long long iter = 0;

    auto BestSolution(x);
    int bestObjVal = computeSolutionObjective(n, a, x);

    //start algo
    while ((clock() < timeToFinish) && iter < grasp_iter_limit) {

        auto grasp_iter_start_time = clock();

        int currentObjVal;
        constructionHeuristic(n, m, x, currentObjVal, alpha, a);
        bool improve;
        localSearch(n, m, x, currentObjVal, a, improve, timeToFinish, jumplevel);

        if (!((iter + 1) % 5000000))
            std::cout << "\n" << std::setw(16) << "Method"  << std::setw(10) << "Iter" << std::setw(13) << "Time" << std::setw(12) <<
                      "ObjVal" << std::setw(20) << "Improves on best?" << std::setw(14) << "Best\n" << std::endl;

        if (!((iter + 1) % 1000000) || (currentObjVal > bestObjVal))
            std::cout << std::setw(16) << "GRASP" << std::setw(10) << iter <<
                  std::setw(10) << (int) ((clock() - grasp_iter_start_time) / double(CLOCKS_PER_SEC / 1000) ) <<
                  " ms" << std::setw(12) << currentObjVal << std::setw(20) << ((currentObjVal > bestObjVal) ? "YES" : "NO") <<
                  std::setw(14) << bestObjVal << std::endl;


        if (currentObjVal > bestObjVal) {
            BestSolution = x;
            bestObjVal = currentObjVal;
        }
        iter++;
    }


    x = BestSolution;

    return iter;
}



int main(int argc, char **argv) {


    std::vector<std::vector<int>> a;
    int n, m;

    int grasp_time_limit = GRASP_TIME_LIMIT;
    long long grasp_iter_limit = GRASP_ITER_LIMIT;
    double alpha = ALPHA;
    int jump_level = JUMP_LEVEL;

#ifdef LOG_TO_FILE
    freopen("log_Improvements.out","a+", stdout);
    printf("\n\n\n\n");         // spacing between different logs

    for (int i = 0; i < argc; ++i) {
        printf(argv[i]);
        printf("  ");
    }
    printf("\n\n");
#endif

    if (argc > 1)   {
        for (int i = 1; i < argc - 1; ++i)


            if (!strcmp(argv[i], "-grasptime")) {
                grasp_time_limit = atoi(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-graspiter")) {
                grasp_iter_limit = atol(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-alpha")) {
                alpha = atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-jumplevel")) {
                jump_level = atoi(argv[i + 1]);
                if (jump_level < 0 || jump_level > 2) {
                    std::cerr << "Wrong Jump Level Specification!";
                    exit(444);
                }

            }
    }


    std::cout << std::setprecision(4) << std::fixed;


    int best_obj;

    for (auto filename: FILES) {
        // initialize
        if (readInstance(filename.insert(0, "../data/"), a, n, m)) {

            Solution x(n); //solution, -1 by default means not processed by the algorithm
            int objVal;
            double time_taken_to_construct = constructionHeuristic(n, m, x, objVal, 1., a);


            Solution initial_solution(x);

            double localSearchTime;
            int init_obj = best_obj = computeSolutionObjective(n, a, initial_solution);
            assert(init_obj == objVal);


            std::cout << std::setw(16) << "Method" << std::setw(10) << "Iter" << std::setw(13) << "Time" << std::setw(12) <<
                      "ObjVal" << std::setw(20) << "Improves on best?" << std::setw(14) << "Best" << std::endl;

            std::cout << std::setw(16) << "Construct" << std::setw(10) << "" <<
                      std::setw(10) << (int) time_taken_to_construct <<
                      " ms" << std::setw(12) << init_obj << std::setw(20) << "N/A" <<
                      std::setw(14) << best_obj << std::endl;

            int graspIterations = GRASP(n, m, x, a, alpha, jump_level, grasp_time_limit, grasp_iter_limit);


            if (graspIterations >= grasp_iter_limit)
                printf("\n\nIteration limit reached.\n");
            else
                printf("\n\nTime limit reached.\n");

            printf("Best Objective Value is %d\n", computeSolutionObjective(n, a, x));
            printSolution(n, x);

        }
    }

    return 0;
}


