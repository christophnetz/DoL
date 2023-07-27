#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <numeric>
#include <fstream>
#include <Eigen>

#include "rndutils.hpp"

using namespace std;
using namespace Eigen;

std::mt19937_64 rng;
int Tmax = 120000;
double lambda = 0.9; //forgetting, lim-> 10
double updaterate = 0.2;
double mort = 0.006; //0.006
double birthrate = 0.08; //0.2 //6, 0.2, 0.06
double phi = 0.0; // interaction term for suitability x experience
double f1 = 0.5;
double f2 = 0.0;


Matrix<double, 1, 3> suitability(const double size) {

  Matrix<double, 1, 3> suit = {
    //1.0 / (1.0 + exp((size - 15) / 4.0)),
    //0.8 / ((1.0 + exp((20 - size) / 3.0))),
    //3.0 / ((1.0 + exp((35 - size) / 4.0)) * (1.0 + exp((size - 35) / 4.0)))
    //    1.0 / (1.0 + exp((25 - size) / 4.0))

        1.0 / (1.0 + exp((size - 15) / 4.0)),
    3.0 / ((1.0 + exp((29 - size) / 4.0)) * (1.0 + exp((size - 29) / 4.0))),
    3.0 / ((1.0 + exp((33 - size) / 4.0)) * (1.0 + exp((size - 33) / 4.0)))
  };


  return suit;
}




double sd(Eigen::ArrayXd vec) {
  return std::sqrt((vec - vec.mean()).square().sum() / (vec.size()));
}




struct ind
{
  ind() : size(1.0), ID(1), updates(0), itask(0) {
    m_experience = { 0.0, 0.0, 0.0 };
    m_task = { 0.33, 0.33, 0.33 };
  };


  ind(int i) : size(1.0), updates(0), itask(0), ID(i) {
    m_experience = { 0.0, 0.0, 0.0 };
    m_task = { 0.33, 0.33, 0.33 };

  };

  double size;
  int updates;
  int ID;
  Matrix<double, 1, 3> m_experience;
  Matrix<double, 1, 3> m_task;
  int itask;

  void experience() {

    for (int i = 0; i < m_experience.size(); ++i) {
      if (i == itask) {
        m_experience(0, itask) += (0.01 * (10.0 - (10.0 * m_experience[itask])));

      }
      else {
        m_experience *= lambda;
      }
    }


  }

  void task(Matrix<double, 1, 3> labour, Matrix<double, 1, 3> count) {


    // mismatch between labour components;
    // sensitivity to mismatch
    Eigen::Index   minIndex;
    if (count.minCoeff(&minIndex) == 0) {
      itask = minIndex;
    }
    else {
      m_task = suitability(size) + m_experience + phi * (suitability(size).array() * m_experience.array()).matrix();

      // one of two decision mechanics: 1) choose task of greatest deficit, or 2) choose task of greatest efficiency gain


      Matrix<double, 1, 3> avgeffectdiff = (m_task.array() / m_task.sum() + f1 * (1.0 - labour.array() / labour.sum()) + f2 * (1.0 - labour.array() / labour.sum()) * m_task.array() / m_task.sum()).matrix();
        
        
      Eigen::Index   maxIndex;
      avgeffectdiff.colwise().sum().maxCoeff(&maxIndex);
      itask = maxIndex;

    }

    updates++;
  }

  void labour(Matrix<double, 1, 3>& labour, Matrix<double, 1, 3>& count) {

    m_task = suitability(size) + m_experience + phi * (suitability(size).array() * m_experience.array()).matrix();


    Matrix<double, 1, 3> mask = { 0.0, 0.0, 0.0 };
    mask(0, itask) = 1.0;

    labour += (mask.array() * (m_task).array()).matrix();//

    count(0, itask) += 1.0;

  }
};

double avg_size(std::vector<ind> v) {
  double totalsize = 0;
  for (int i = 0; i < v.size(); i++) {
    totalsize += v[i].size;
  }
  return totalsize / static_cast<double>(v.size());
}

void deaths(vector<ind>& pop) {

  std::binomial_distribution<int> d_dist(pop.size(), mort);
  int deaths = d_dist(rng);
  if (pop.size() - deaths < 5) {  //setting a minimal population size
    deaths = pop.size() - 5;
  }
  if (deaths > 0) {

    std::shuffle(std::begin(pop), std::end(pop), rng);
    pop.resize(pop.size() - deaths);
  }
}

void births(vector<ind>& pop, int& ID, const Matrix<double, 1, 3> labour, const Matrix<double, 1, 3> count) {

  int offspring = std::poisson_distribution<int>(birthrate)(rng);
  double popsize = static_cast<double> (pop.size());
  for (int i = 0; i < offspring; ++i) {
    pop.push_back(ind(ID));
    pop.back().task(labour, count);
    ID = (ID + 1);
  }
}



int main() {

  std::ofstream ofs1("output1small.txt", std::ofstream::out);
  std::ofstream ofs2("output2small.txt", std::ofstream::out);
  ofs1 << "t\tpopsize" << "\t" << "avgbodysize" << "\t" << "sdneed\teggcare\tdigging\tdefense" << std::endl;
  ofs2 << "t\tID\tsize\tcurrent_task\tupdates\texp_egg\texp_digg\texp_def\ttask_egg\ttask_digg\ttask_def" << std::endl;

  Matrix<double, 1, 3> labour = { 3.0, 3.0, 3.0 };
  Matrix<double, 1, 3> counts = { 10.0, 10.0, 10.0 };

  vector<ind> pop(30);
  int ID = 0;

  rng.seed(5555);

  for (int i = 0; i < pop.size(); ++i) {
    pop[i].size = static_cast<double>(i) + 1.0;
  }

  int next_t = 1;
  // discrete times, but don't update simultaneously, instead draw individuals at random one after the other
  for (double t = 0.0; t < Tmax; ) {

    double delta_t = std::exponential_distribution<double>(static_cast<double>(pop.size()) * updaterate)(rng);

    while (delta_t + t > static_cast<int>(next_t)) {

      births(pop, ID, labour, counts);
      deaths(pop);

      labour = { 0.0, 0.0, 0.0 }; //reset labour vector
      counts = { 0.0, 0.0, 0.0 }; //reset counts 

      for (int i = 0; i < pop.size(); i++) {
        pop[i].labour(labour, counts);
        pop[i].experience();
        pop[i].size *= (1.0 + 0.2 * exp(-0.15 * pop[i].size));
      }

      // Output
      ofs1 << next_t << "\t" << pop.size() << "\t" << avg_size(pop) << "\t" << sd(labour.array()) << "\t" << labour << std::endl;

      if (next_t > 99000.0) { //90000
        for (int i = 0; i < pop.size(); ++i) {

          ofs2 << next_t << "\t" << pop[i].ID << "\t" << pop[i].size << "\t" << pop[i].itask << "\t" << pop[i].updates << "\t" << pop[i].m_experience << "\t" << pop[i].m_task << std::endl;

        }
      }
      next_t++;
    }



    int sel_i = std::uniform_int_distribution<int>(0, pop.size() - 1) (rng);

    pop[sel_i].task(labour, counts);

    t += delta_t;

  }

  return 0;
}


