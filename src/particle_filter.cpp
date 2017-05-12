/*
 * particle_filter.cpp
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // Set the number of particles
    num_particles = 200;
    
    // Set state noise Gaussians
    default_random_engine gen;
    normal_distribution<double> N_x_init(0, std[0]);
    normal_distribution<double> N_y_init(0, std[1]);
    normal_distribution<double> N_theta_init(0, std[2]);
    double n_x, n_y, n_theta;
    
    // Iterate particle generation
    for (int i=0; i<num_particles; ++i) {
        
        // Generate noise
        n_x = N_x_init(gen);
        n_y = N_y_init(gen);
        n_theta = N_theta_init(gen);
        
        // Generate a particle
        Particle p;
        
        // Set the particle's parameters and add noise
        p.id = i;
        p.x = x + n_x;
        p.y = y + n_y;
        p.theta = theta + n_theta;
        p.weight = 1;
        
        // Append this particle to the set
        particles.push_back(p);
    }

    // Set the flag
    is_initialized = true;
    
    // Append the weight to weights vector
    fill(weights.begin(), weights.end(), 1);
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    // Set state noise Gaussian
    default_random_engine gen;
    normal_distribution<double> N_x(0, std_pos[0]);
    normal_distribution<double> N_y(0, std_pos[1]);
    normal_distribution<double> N_theta(0, std_pos[2]);
    double n_x, n_y, n_theta;
    
    // Predict for each particle
    for (int i=0; i<num_particles; ++i) {
        
        // Current state
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        
        // Predict with given velocity, yaw rate and delta time
        if (fabs(yaw_rate) > 0.0001) {    // Nonlinear motion
            double theta_f = fmod(theta + yaw_rate * delta_t, 2.0 * M_PI);
            x += velocity / yaw_rate * (sin(theta_f) - sin(theta));
            y += velocity / yaw_rate * (cos(theta) - cos(theta_f));
            theta = theta_f;
        } else {    // Straight motion
            x += velocity * cos(theta) * delta_t;
            y += velocity * sin(theta) * delta_t;
        }
        
        // Generate noise
        n_x = N_x(gen);
        n_y = N_y(gen);
        n_theta = N_theta(gen);
        
        // Add noise and set the predicted state
        particles[i].x = x + n_x;
        particles[i].y = y + n_y;
        particles[i].theta = theta + n_theta;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
 
    // Find the nearest landmark for all observations
    for (int i=0; i<observations.size(); ++i) {
        
        // Initialize minimum distance from observation to landmarks
        double dist_min = 0.0;
        // Id of matched landmark
        int id_m = 0;
        
        // Iterate over all landmarks
        for (int j=0; j<predicted.size(); ++j) {
            
            // Find the distance from observation to landmark
            double dist_temp = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            
            // Set initial minimum distance and id of matched landmark
            if (j == 0) {
                dist_min = dist_temp;
                id_m = predicted[j].id;
            }
            
            // Find the minimum distance and corresponding landmark id
            if (dist_temp < dist_min) {
                dist_min = dist_temp;
                id_m = predicted[j].id;
            }
        }
        
        // Set the matched landmark id to the observation
        observations[i].id = id_m;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33.)
	//   http://planning.cs.uiuc.edu/node99.html
    //   Wiki: rotation matrix
    //   https://en.wikipedia.org/wiki/Rotation_matrix
    
    // Declare a temporary weights vector
    vector<double> weights_temp;
    
    // Iterate over all particles
    for (int n=0; n<num_particles; ++n) {

        /*******************************************************
                Transform observations to map coordinates
         *******************************************************/
        
        // Set values for better readability
        const double x = particles[n].x;
        const double y = particles[n].y;
        const double theta = particles[n].theta;
        
        // Declare a set of transformed observations
        vector<LandmarkObs> observations_trans;
        
        // Declare transformed values
        double x_trans;
        double y_trans;
        
        // Transform all observations
        for (int i=0; i<observations.size(); ++i) {
            
            // Transform: rotation and translation
            x_trans = x + observations[i].x*cos(theta) - observations[i].y*sin(theta);
            y_trans = y + observations[i].x*sin(theta) + observations[i].y*cos(theta);
            
            // Add the transformed observation to the set
            LandmarkObs obs;
            obs.x = x_trans;
            obs.y = y_trans;
            observations_trans.push_back(obs);
        }
        
        /*******************************************************
                Find map landmarks within sensor range
         *******************************************************/

        // Declare a set of landmarks within sensor range
        vector<LandmarkObs> predicted;
        
        // Declare landmark coordinates and id
        double x2, y2;
        int id;
        
        // Iterate over all landmarks in map
        for (int m=0; m<map_landmarks.landmark_list.size(); ++m) {
            
            // Find the coordinates and id of the landmark
            x2 = map_landmarks.landmark_list[m].x_f;
            y2 = map_landmarks.landmark_list[m].y_f;
            id = map_landmarks.landmark_list[m].id_i;
            
            // Check if the landmark is within sensor range
            if (dist(x, y, x2, y2) <= sensor_range) {
                
                // Add the within-range landmark to the set
                LandmarkObs landmark;
                landmark.id = id;
                landmark.x = x2;
                landmark.y = y2;
                predicted.push_back(landmark);
            }
        }
        
        /*******************************************************
                Match observations with landmarks
         *******************************************************/
        
        dataAssociation(predicted, observations_trans);
        
        /*******************************************************
                Update importance weight
         *******************************************************/
        
        // Declare variables for better readability
        double x_obs, y_obs, x_landmark, y_landmark;
        int id_landmark;
        
        // Declare importance weight
        double prob = 1.0;
        
        // Iterate over all observations
        for (int i=0; i<observations_trans.size(); ++i) {
            
            // Get the coordinates of corresponding observation and landmark
            x_obs = observations_trans[i].x;
            y_obs = observations_trans[i].y;
            id_landmark = observations_trans[i].id - 1; // landmap list indicing from 0
            x_landmark = map_landmarks.landmark_list[id_landmark].x_f;
            y_landmark = map_landmarks.landmark_list[id_landmark].y_f;
            
            // Calculate the importance weight (measurement probability)
            // using normdf method defined in helper_functions.h
            prob *= normpdf(x_obs, x_landmark, std_landmark[0]);
            prob *= normpdf(y_obs, y_landmark, std_landmark[1]);
        }
        
        // Update weight
        particles[n].weight = prob;
        
        // Add the weight to the temporary weights vector
        weights_temp.push_back(prob);
    }
    
    // Update the weights vector of particle filter
    weights = weights_temp;

    /*******************************************************
            Normalize importance weights
     *******************************************************/

    // In fact this step can be saved, since discrete_distribution
    // in the resample method will draw sample as per normalized weights.
    //weights = normalize_vector(weights);
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    // Define a distribution in which sample indices, a series of integers,
    // will be drawn as per the weights vector.
    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<int> d(weights.data(), weights.data()+weights.size());
    
    // Declare resampled particles
    vector<Particle> particles_resampled;
    
    // Resample
    for (int n=0; n<num_particles; ++n) {
        particles_resampled.push_back(particles[d(gen)]);
    }
    particles = particles_resampled;
    
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
