#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

namespace ublas = boost::numeric::ublas;

namespace constants{
const double sigma = 0.86;
const double epsilon = 0.0001;
}

using ublas_vector = ublas::vector<double, std::vector<double>>;

class Vector3D : public ublas_vector{
 public:
  using ublas_vector::ublas_vector;
  Vector3D& operator=(std::initializer_list<double> list){
    std::move(list.begin(), list.end(), this->begin());
    return *this;
  }
};

struct Particle{
  Particle() : position(3), position_old(3),
               velocity(3), velocity_old(3),
               force(3), force_old(3){}
  Vector3D position;
  Vector3D position_old;
  Vector3D velocity;
  Vector3D velocity_old;
  Vector3D force;
  Vector3D force_old;
};

class Simulation {
 public:
  void Run() {
    InitializeParticles();

    Save("../data/data_0.csv"); // TODO: Create directory by program
    for (int i{1}; i <= steps_number; ++i) {
      std::cout << "Step: " << i << '\n';
      Step();
      Save("../data/data_" + std::to_string(i) + ".csv");

    }

  }

  void Save(const std::string &path) {
    std::ofstream file(path, std::ios::out);
    for (int i{}; i < particles_number; ++i) {
      file << particles[i].position[0] << ','
           << particles[i].position[1] << ','
           << particles[i].position[2] << ','
           << particles[i].velocity[0] << ','
           << particles[i].velocity[1] << ','
           << particles[i].velocity[2] << '\n';
    }
  }

 private:
  const std::vector<double> box_size{25, 25, 25};
  const int particles_number = 100;
  const int steps_number = 100000;
  const double dt = 0.001;

  // max distance when force/energy computed yet
  const double cutoff = std::min(5.0*constants::sigma,
                                 *std::min_element(box_size.cbegin(), box_size.cend()));

  std::vector<Particle> particles;

  void InitializeParticles() {
    particles.resize(particles_number);
    int particles_count = 0;

    auto num_per_length = pow(particles_number, 1.0 / 3) + 1;

    std::uniform_real_distribution<double> distribution(-10, 10);
    std::random_device rd;
    std::mt19937 engine(rd());

    auto dx = box_size[0] / num_per_length;
    auto dy = box_size[1] / num_per_length;
    auto dz = box_size[2] / num_per_length;
    for (int i0{1}; i0 <= num_per_length; ++i0) {
      for (int i1{1}; i1 <= num_per_length; ++i1) {
        for (int i2{1}; i2 <= num_per_length; ++i2) {
          if (particles_count < particles_number) {
            auto &p = particles[particles_count];
            p.position_old = p.position = {dx * i0, dy * i1, dz * i2};

            p.velocity_old = p.velocity = {distribution(engine), distribution(engine), distribution(engine)};
            particles_count++;
          }
        }
      }
    }
  }
  void Step() {
    UpdateForces();
    UpdateVelocities();
    UpdatePositions();
    PeriodicCondition();

  }

  double potential_energy(int i, int j) {
    auto dr = BoundaryCondition(i, j);
    auto dist2 = ublas::inner_prod(dr, dr);
    if (dist2 > cutoff*cutoff){
      return 0.0;
    }
    auto r_inv2 = constants::sigma*constants::sigma/dist2;
    auto r_inv6 = r_inv2*r_inv2*r_inv2;
    auto r_inv12 = r_inv6*r_inv6;

    return 4 * constants::epsilon * (r_inv12 - r_inv6);
  }

  Vector3D force(int i, int j) {
//    auto dr = particles[i].position - particles[j].position;
    auto dr = BoundaryCondition(i, j);
    auto dist2 = ublas::inner_prod(dr, dr);
    if (dist2 > cutoff*cutoff){
      return Vector3D({0.0, 0.0, 0.0});
    }
    auto r_inv2 = constants::sigma*constants::sigma/dist2;
    auto r_inv6 = r_inv2*r_inv2*r_inv2;
    auto r_inv12 = r_inv6*r_inv6;

    return 24 * constants::epsilon*r_inv2 * (2 * r_inv12 - r_inv6) * dr;
  }

  void UpdateForces() {
    for (int i{}; i < particles_number; ++i) {
      particles[i].force_old = particles[i].force;
      particles[i].force = {0.0, 0.0, 0.0};
    }

    for (int i{}; i < particles_number; ++i) {
      for (int j{i + 1}; j < particles_number; ++j) {
        auto f = force(i, j);
        particles[i].force += f;
        particles[j].force -= f;
      }
    }
  }

  void UpdateVelocities() {
    for (int i{}; i < particles_number; ++i) {
      particles[i].velocity_old = particles[i].velocity;
      particles[i].velocity += 0.5 * (particles[i].force + particles[i].force_old) * dt;
    }
  }

  void UpdatePositions() {
    for (int i{}; i < particles_number; ++i) {
      particles[i].position_old = particles[i].position;
      particles[i].position +=
          particles[i].velocity * dt + 0.5 * particles[i].force * dt * dt;
    }
  }
  void PeriodicCondition() {
    // if particles leaves central cell it will be returned from the opposite side
    for (int i{}; i < particles_number; ++i) {
      auto &p = particles[i];
      for (int j{}; j<3; ++j) {
        if (p.position[j] < 0) {
          p.position[j] += box_size[j];
        } else if (p.position[j] >= box_size[j]) {
          p.position[j] -= box_size[j];
        }
      }
    }
  }

  Vector3D BoundaryCondition(int i, int j){
    // particle 'i' interacts only with nearest particle 'j' or its clone (along each axis)
    Vector3D dr(particles[i].position - particles[j].position);

    for (int k{}; k<3; ++k){
      if (dr[k] < -0.5*box_size[k]){
        dr[k] += box_size[k];
      }
      if (dr[k] > 0.5*box_size[k]){
        dr[k] -= box_size[k];
      }
    }

    return dr;
  }

};

int main(int argv, char** argc){
  Simulation simulation;

  simulation.Run();

  return EXIT_SUCCESS;
}