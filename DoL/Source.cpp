#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <numeric>
#include <fstream>
#include <iostream>
#include <Eigen>
#include "rndutils.hpp"
#include "parameters.h"
#include "cmd_line.h"

using namespace std;
using namespace Eigen;

std::mt19937_64 rng;



Matrix<double, 1, 2> baseline() {

  Matrix<double, 1, 2> base = {
    1.0,
    1.0
  };

  return base;
}


double sd(Eigen::ArrayXd vec) {
  return std::sqrt((vec - vec.mean()).square().sum() / (vec.size()));
}


struct ind
{
  ind() : size(1.0), ID(1), updates(0), itask(0), task_changes(0) {
    m_experience = { 0.0, 0.0 };
    m_task = { 0.00, 0.00 };
  };


  ind(int i) : size(1.0), updates(0), itask(0), task_changes(0), ID(i) {
    m_experience = { 0.0, 0.0 };
    m_task = { 0.00, 0.00 };

  };

  double size;
  int updates;
  int ID;
  Matrix<double, 1, 2> m_experience;
  Matrix<double, 1, 2> m_task;
  int itask;
  int task_changes; 

  void experience(double forget_rate, double learning_rate) {

    for (int i = 0; i < m_experience.size(); ++i) {
      if (i == itask) {
        m_experience(0, i) += (learning_rate * (1.0 - (1.0 * m_experience[itask])));
      }
      else {
        m_experience(0, i) *= (1.0 - forget_rate);
      }
    }
  }


  void task(Matrix<double, 1, 2>& labour, Matrix<double, 1, 2> count, param_t params) {

    labour[itask] -= m_task[itask];

    m_task = baseline() + m_experience;

    Matrix<double, 1, 2> avgeffectdiff = ((1.0 - params.f1) * m_task.array() / m_task.sum() + params.f1 * (1.0 - labour.array() / labour.sum())).matrix();


    Eigen::Index   maxIndex;
    avgeffectdiff.colwise().sum().maxCoeff(&maxIndex);
    if (itask != maxIndex) { 
      itask = maxIndex; 
      task_changes++; 
    }
    itask = maxIndex;
    labour[itask] += m_task[itask];

    updates++;
  }


  void labour(Matrix<double, 1, 2>& labour, Matrix<double, 1, 2>& count, param_t params) {

    m_task = baseline() + m_experience;

    Matrix<double, 1, 2> mask = { 0.0, 0.0 };
    mask(0, itask) = 1.0;

    labour += (mask.array() * (m_task).array()).matrix();

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

void deaths(vector<ind>& pop, double mortrate, std::ofstream& ofs) {

  std::binomial_distribution<int> d_dist(pop.size(), mortrate);
  int deaths = d_dist(rng);
  if (pop.size() - deaths < 2) {  //setting a minimal population size

    deaths = pop.size() - 1;
  }

  if (deaths > 0) {

    std::shuffle(std::begin(pop), std::end(pop), rng);

    for (int i = 0; i < deaths; i++) {  
      ofs << pop.back().task_changes <<" ";
      pop.pop_back();
    } 
  }
}

void births(vector<ind>& pop, int& ID, bool idactive, Matrix<double, 1, 2>& labour, const Matrix<double, 1, 2> count, param_t params) {

  int offspring = std::poisson_distribution<int>(params.birthrate)(rng);
  double popsize = static_cast<double> (pop.size());
  for (int i = 0; i < offspring; ++i) {
    pop.push_back(ind(ID));

    pop.back().task(labour, count, params);

    if (idactive) {
      ID = (ID + 1);
    }
  }
}



void run_sim(param_t params) {

  int _N = 0;
  double a;
  double b;
  double avgsd = 0.0;
  double avgperf = 0.0;


  string str = params.outdir + "_summary.txt";
  string str1 = params.outdir + "_ind.txt";
  string str2 = params.outdir + "_changes.txt";

  std::ofstream ofssum(str, std::ofstream::out | std::ofstream::app);
  std::ofstream ofs2(str1, std::ofstream::out | std::ofstream::app);
  std::ofstream ofs_changes(str2, std::ofstream::out | std::ofstream::app);
  ofs_changes << "f1 "<< params.f1<< "\tbirthrate\t" << params.birthrate << "\learning\t" << params.learning_rate  << " "; 
  //ofssum << "samples\tavgperf\tavgsd\tf1\tbirthrate\tlearn\tforget\t" << std::endl;
  ofs2 << "f1\tt\tID\tcurrent_task\texp_egg\texp_digg\texp_def\tnr_changes" << std::endl;


  Matrix<double, 1, 2> labour = { 3.0, 3.0 };
  Matrix<double, 1, 2> counts = { 10.0, 10.0 };

  vector<ind> pop(30);
  int ID = 0;
  int ind_counter = 0;
  if (params.seed == 0) {
    static std::random_device rd{}; 
    auto seed = rd();
    rng.seed(seed); // seed the engine
    cout << "Seed: "<< seed << endl;
  }
  else {
    rng.seed(params.seed); 
  }

  for (int i = 0; i < pop.size(); ++i) {
    pop[i].size = static_cast<double>(i) + 1.0;
  }

  int next_t = 1;
  // discrete times, but don't update simultaneously, instead draw individuals at random one after the other
  for (double t = 0.0; ID < params.Tmax; ) {

    double delta_t = std::exponential_distribution<double>(static_cast<double>(pop.size()) * params.updaterate)(rng);

    while (delta_t + t > next_t) {

      births(pop, ID, next_t > 5000, labour, counts, params);
      deaths(pop, params.mort, ofs_changes);

      labour = { 0.0, 0.0 }; //reset labour vector
      counts = { 0.0, 0.0 }; //reset counts 

      for (int i = 0; i < pop.size(); i++) {
        pop[i].size *= (1.0 + 0.2 * exp(-0.15 * pop[i].size));
        pop[i].experience(params.forget_rate, params.learning_rate);
        pop[i].labour(labour, counts, params);
      }


      if (next_t > 5000.0) { 

        _N++;
        a = 1.0 / static_cast<double>(_N);
        b = 1.0 - a;
        avgperf = a * labour.sum()/static_cast<double>(pop.size()) + b * avgperf;
        avgsd = a * sd(labour.array()/ labour.sum()) + b * avgsd;


        for (int i = 0; i < pop.size(); ++i) {
          if (pop[i].ID > 0) {
            ofs2 << params.f1 << "\t" << next_t << "\t" << pop[i].ID << "\t" << pop[i].itask << "\t" << pop[i].m_experience << "\t" << pop[i].task_changes << std::endl;
          }
        }
      }
      next_t++;
    }



    int sel_i = std::uniform_int_distribution<int>(0, pop.size() - 1) (rng);

    pop[sel_i].task(labour, counts, params);

    t += delta_t;

  }

  ofssum << _N << "\t" << avgperf << "\t" << avgsd << "\t" << params.f1 << "\t" << params.birthrate << "\t"
    << params.learning_rate << "\t" << params.forget_rate<< endl;

  ofssum.close();
  ofs_changes << std::endl;
  ofs_changes.close();
  ofs2.close();
}

int main(int argc, const char** argv) {

  try {
    auto param = param_t{};
    cmd::cmd_line_parser clp(argc, argv);
    auto config = clp.optional_val("config", std::string{});
    if (!config.empty()) clp.append(config_file_parser(config));
    clp.optional("Tmax", param.Tmax);
    clp.optional("forget_rate", param.forget_rate);
    clp.optional("learning_rate", param.learning_rate);
    clp.optional("updaterate", param.updaterate);
    clp.optional("mort", param.mort);
    clp.optional("birthrate", param.birthrate);
    clp.optional("f1", param.f1);
    clp.optional("seed", param.seed);
    clp.optional("outdir", param.outdir);


    auto unknown = clp.unrecognized();
    if (!unknown.empty()) {
      for (const auto& arg : unknown) {
        std::cerr << "unknown argument '" << arg << "'\n";
      }
      return 1;
    }

    run_sim(param);

    return 0;
  }

  catch (cmd::parse_error& err) {
    std::cerr << "\nParameter trouble: " << err.what() << '\n';
  }
  catch (std::exception& err) {
    std::cerr << "\nExeption caught: " << err.what() << '\n';
  }
  catch (...) {
    std::cerr << "\nUnknown exeption caught\n";
  }

  return 0;
}
