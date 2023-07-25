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
int Tmax = 100000;
double lambda = 0.9; //forgetting, lim-> 10
double updaterate = 0.2;
double mort = 0.006; //0.006
double birthrate = 0.06; //0.2 //6, 0.2, 0.06
double phi = 0.0; // interaction term for suitability x experience


Matrix<double, 1, 3> suitability(const double size) {

  Matrix<double, 1, 3> suit = { 
    1.0 / (1.0 + exp((size - 15) / 4.0)), 
    3.0 / ((1.0 + exp((25 - size) / 4.0)) * (1.0 + exp((size - 25) / 4.0))),
    1.0 / (1.0 + exp((25 - size) / 4.0)) };
  
 
  return suit;
}




double sd(Eigen::ArrayXd vec) {
  return std::sqrt((vec - vec.mean()).square().sum() / (vec.size()));
}




struct ind
{
  ind() : size(1.0), ID(1), updates(0) { 
    m_experience = { 0.0, 0.0, 0.0 };
    m_task = { 0.33, 0.33, 0.33 };
  };


  ind(int i) : size(1.0), updates(0), ID(i) { 
    m_experience = { 0.0, 0.0, 0.0 }; 
    m_task = { 0.33, 0.33, 0.33 };

  };

  double size;
  int updates;
  int ID;
  Matrix<double, 1, 3> m_experience;
  Matrix<double, 1, 3> m_task;

  void experience() {

    m_experience = m_task + lambda * m_experience;
  }

  void task( Matrix<double, 1, 3> labour, double popsize) {

    /*
    get min of labour - pour in individual labour until equalized, distribute rest evenly
    how much of their 'budget' follows intrinsic strength, how much follows balance requirement?

    some scale may be needed, population labour vs individual contribution
    skill = suitability(size, i) + m_experience[i]

         

    double time = 0.0;
    2nd labourdist - min(labourdist) = x;
    max(labourdist) - 2nd labourdist = y;

    if (suitability(size, i) + m_experience[i]) < x
      task[min] = 1, others 0;
    else
      task[min] = x / (suitability(size, i) + m_experience[i]);

    if (suitability(size, i) + m_experience[i]) < 2 * y{
      task[min] += ;
    task[2nd] += ;
    }
    else {
      task[min] += x / (suitability(size, i) + m_experience[i]);
      task[2nd] += x / (suitability(size, i) + m_experience[i]);
    }

    for (int i = 0; i < m_task.size(); i++) {
      task[i] += (suitability(size, i) + m_experience[i]);


    }

      
 */


    m_task =  suitability(size) + m_experience/10.0 - labour/popsize + phi * (suitability(size).array() * m_experience.array() / 10.0).matrix() ; //
    m_task.array() += 1;  //really this is (N - labour), but due to eigen library better written this way
    updates++;
    m_task = m_task / m_task.sum();
  }

  void labour(Matrix<double, 1, 3>& labour) {

    labour += (m_task.array() * (suitability(size) + m_experience/10.0  + phi * (suitability(size).array() * m_experience.array() / 10.0).matrix()  ).array()).matrix();//
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

void births(vector<ind>& pop, int& ID, const Matrix<double, 1, 3> labour) {
  
  int offspring = std::poisson_distribution<int> (birthrate)(rng);
  double popsize = static_cast<double> (pop.size());
  for (int i = 0; i < offspring; ++i) {
    pop.push_back(ind(ID));
    pop.back().task(labour, popsize);
    ID = (ID + 1);
  }
}




int main() {

  std::ofstream ofs1("output1.txt", std::ofstream::out);
  std::ofstream ofs2("output2.txt", std::ofstream::out);
  ofs1 << "t\tpopsize" << "\t" << "avgbodysize" << "\t" << "sdneed\teggcare\tdigging\tdefense" << std::endl;
  ofs2 << "t\tID\tsize\tupdates\texp_egg\texp_digg\texp_def\ttask_egg\ttask_digg\ttask_def" << std::endl;

  Matrix<double, 1, 3> labour = { 0.0, 0.0, 0.0 };

  vector<ind> pop(30);
  int ID = 0;

  for (int i = 0; i < pop.size(); ++i) {
    pop[i].size = static_cast<double>(i) + 1.0;
  }

  int next_t = 1;
  // discrete times, but don't update simultaneously, instead draw individuals at random one after the other
  for (double t = 0.0; t < Tmax; ) {

    double delta_t = std::exponential_distribution<double>(static_cast<double>(pop.size()) * updaterate)(rng);

    while (delta_t + t > static_cast<int>(next_t)) {
      
      births(pop, ID, labour);
      deaths(pop);

      labour = { 0.0, 0.0, 0.0 }; //reset labour vector

      for (int i = 0; i < pop.size(); i++) {
        pop[i].labour(labour);
        pop[i].experience();
        pop[i].size *= (1.0 + 0.2 * exp(-0.15 * pop[i].size));
      }

      // Output
      ofs1 << next_t << "\t" << pop.size() << "\t" << avg_size(pop) << "\t" << sd(labour.array()) << "\t" << labour << std::endl;

      if (next_t > 90000.0 && next_t % 10 == 0) { //90000
        for (int i = 0; i < pop.size(); ++i) {

          ofs2 << next_t << "\t" << pop[i].ID << "\t" << pop[i].size << "\t" << pop[i].updates << "\t" << pop[i].m_experience << "\t" << pop[i].m_task << std::endl;

        }
      }
      next_t++;
    }


    int sel_i = std::uniform_int_distribution<int>(0, pop.size() - 1) (rng);

    pop[sel_i].task(labour, static_cast<double>(pop.size()));

    t += delta_t;

  }

  return 0;
}


