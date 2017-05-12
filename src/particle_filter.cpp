/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
bool debug = false;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	// x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if(debug)
		std::cout << "------------------------INITIALIZATION------------------------" << std::endl;
	num_particles = 200;

	std::default_random_engine gen;
	std::normal_distribution<double> distribution_x(x, std[0]);
	std::normal_distribution<double> distribution_y(y, std[1]);
	std::normal_distribution<double> distribution_t(theta, std[2]);
	for (int i = 0; i < num_particles; i++)
	{
		Particle p =
		{
			i,					   //particle ID
			distribution_x(gen),   //position x
			distribution_y(gen),   //position y
			distribution_t(gen),   //theta
			1.0            		   //weights
		};
		particles.push_back(p);
		weights.push_back(p.weight);
		if(debug)
			std::cout << "Init Particle: "<< particles[i].id << "," << particles[i].x << "," << particles[i].y << "," << particles[i].theta << "," << particles[i].weight << std::endl;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//Generate normal distribution to sample random values for noise
	std::default_random_engine gen;
	for (int i=0;i < num_particles; i++)
	{
		if(debug)
		{
			std::cout << "------------------------PREDICTION------------------------" << std::endl;
			std::cout << "Start Prediction: "<< particles[i].x<< "," << particles[i].y << "," << particles[i].theta << std::endl;
		}
		// Check for divide by zero on yaw_rate and switch between motion models
		double theta_delt          = particles[i].theta + (yaw_rate*delta_t);
		if(fabs(yaw_rate) < 0.00001)
		{
			particles[i].x     += velocity * delta_t * cos(particles[i].theta);
			particles[i].y     += velocity * delta_t * sin(particles[i].theta);
			particles[i].theta  = 0;
		}
		else
		{
			particles[i].x     += (velocity / yaw_rate) *  (sin(theta_delt) - sin(particles[i].theta));
			particles[i].y     += (velocity / yaw_rate) *  (cos(particles[i].theta) - cos(theta_delt));
			particles[i].theta += yaw_rate * delta_t;
		}

		// Add noise
		std::normal_distribution<double> distribution_x(0, std_pos[0]);
		std::normal_distribution<double> distribution_y(0, std_pos[1]);
		std::normal_distribution<double> distribution_theta(0, std_pos[2]);

		particles[i].x     = particles[i].x + distribution_x(gen);
		particles[i].y     = particles[i].y + distribution_y(gen);
		particles[i].theta = particles[i].theta + distribution_theta(gen);
		if(debug)
			std::cout << "End Prediction: "<< particles[i].x<< "," << particles[i].y << "," << particles[i].theta << std::endl;
	}//end for loop over all particles

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations)
{
	// Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	// NOT USED - SEE UPDATE WEIGHTS FUNCTION FOR NEAREST NEIGHBOR IMPLEMENTATION
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks)
{
	// Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	//Loop over all particles
	for(int i = 0; i < num_particles; i++)
	{
		if(debug)
			std::cout << "------------------------PROCESSING PARTICLE: "<< i << "------------------------" << std::endl;
		double particle_x     = particles[i].x;
		double particle_y     = particles[i].y;
		double particle_theta = particles[i].theta;

	    // Transform from VEHICLE (observation) to MAP (particle) Coordinates///////////////////////////////////
		if(debug)
			std::cout << "------------------------TRANSFORM OBSERVATIONS:------------------------" << std::endl;
		std::vector<LandmarkObs> vec_transformed_obs;
		for (int j = 0; j < observations.size(); j++)
		{
		  //Rotation and translation between coordinate systems
		  LandmarkObs temp_obs;
		  temp_obs.x  = (observations[j].x * cos(particle_theta)) - (observations[j].y * sin(particle_theta)) + particle_x;
		  temp_obs.y  = (observations[j].x * sin(particle_theta)) + (observations[j].y * cos(particle_theta)) + particle_y;
		  temp_obs.id = -1;
		  vec_transformed_obs.push_back(temp_obs);
		  if(debug)
		  {
			  std::cout << "------------------------" << std::endl;
			  std::cout << "Map Coordinates: "    << temp_obs.x << "," << temp_obs.y << std::endl;
			  std::cout << "Vehicle Coordinates: "<< observations[j].x << "," << observations[j].y << std::endl;
		  }
		}//end loop over transforms

		// Find landmarks within sensor range of the particle to reduce search space////////////////////////////
		if(debug)
			std::cout << "------------------------FIND LANDMARKS IN SENSOR RANGE:------------------------" << std::endl;
		std::vector<LandmarkObs> vec_landmarks_nearby;
		for(int k = 0; k < map_landmarks.landmark_list.size(); k++)
		{
		  double dist_landmark2particle = dist(particle_x, particle_y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
		  //Check if the landmark and particle are within the sensor's range
		  if(dist_landmark2particle <= sensor_range)
		  {
			if(debug)
			{
				std::cout << "Distance Landmark to Particle " << dist_landmark2particle << std::endl;
				std::cout << "Sensor Range " << sensor_range << std::endl;
				std::cout << "Particle: "<< particles[i].x<< "," << particles[i].y << std::endl;
			    std::cout << "Landmark: "<< map_landmarks.landmark_list[k].x_f << "," << map_landmarks.landmark_list[k].y_f << std::endl;
			    std::cout << "Landmark id: "<< map_landmarks.landmark_list[k].id_i << std::endl;
			}
			LandmarkObs temp_landmark;
			temp_landmark.x    = map_landmarks.landmark_list[k].x_f;
			temp_landmark.y    = map_landmarks.landmark_list[k].y_f;
			temp_landmark.id   = map_landmarks.landmark_list[k].id_i;
			vec_landmarks_nearby.push_back(temp_landmark);
		  }
		}//end loop over landmarks

		// Match nearest landmark with transformed observation/////////////////////////////////////////////////
		std::vector<LandmarkObs> nearby_landmark;
		double distance_obs2landm;
		for (int i=0; i<vec_transformed_obs.size(); i++)
		{
			LandmarkObs obs = vec_transformed_obs[i];
			double min_dist = 9999.0;
			int nearby_idx = -1;

			for (int j=0; j<vec_landmarks_nearby.size(); j++)
			{
				distance_obs2landm = dist(obs.x, obs.y, vec_landmarks_nearby[j].x, vec_landmarks_nearby[j].y);
				if (distance_obs2landm < min_dist)
				{
					min_dist      = distance_obs2landm;
					nearby_idx = j;
				}
			}
			if (vec_landmarks_nearby.size() > 0)
			{
				nearby_landmark.push_back(vec_landmarks_nearby[nearby_idx]);
			}
		}//end loop over transformed obs
		if(debug)
		{
			std::cout << "------------------------DATA ASSOCIATION:------------------------" << std::endl;
			for(int h=0; h<vec_transformed_obs.size(); h++)
			{
				std::cout << "nearby_landmark: " << nearby_landmark[h].x << "," << nearby_landmark[h].x << std::endl;
				std::cout << "vec_transformed_obs: " << vec_transformed_obs[h].x << "," << vec_transformed_obs[h].y << std::endl;
				std::cout << "-------------------------------" << std::endl;
			}
		}

		//Compute the weights/////////////////////////////////////////////////////////////////////////////////
		if(debug)
			std::cout << "------------------------COMPUTE WEIGHTS:------------------------" << std::endl;
		double weight = 1.0;

	    for (int j = 0; j < vec_transformed_obs.size(); j++)
	    {
	      LandmarkObs observation = vec_transformed_obs[j];
	      double mu_x  = nearby_landmark[j].x;
	      double mu_y  = nearby_landmark[j].y;
	      if(debug)
	      {
			  std::cout << "-------Match Observation to Prediction-------" << std::endl;
			  std::cout << "Observation: "<< observation.x<< "," << observation.y << std::endl;
			  std::cout << "Prediction: " << mu_x << "," << mu_y << std::endl;
	      }
	      double den      = (2.0*M_PI*std_landmark[0]*std_landmark[1]);
	      double delx     = pow(observation.x-mu_x, 2.0) / (pow(std_landmark[0], 2.0));
	      double dely     = pow(observation.y-mu_y, 2.0) / (pow(std_landmark[1], 2.0));
	      double expterm  = -0.5*(delx + dely);
	      double num      = exp(expterm);
	      double prob_xy  = num/den;
	      weight *= prob_xy;
	      if(debug)
	      {
	    	  std::cout << "-------Weight Computation-------" << std::endl;
			  std::cout<<"delx: "    << delx    << std::endl;
			  std::cout<<"dely: "    << dely    << std::endl;
			  std::cout<<"expterm: " << expterm << std::endl;
			  std::cout<<"num: "     << num     << std::endl;
			  std::cout<<"den: "     << den     << std::endl;
			  std::cout<<"prob_xy: " << prob_xy << std::endl;
			  std::cout<<"weight: "  << weight  << std::endl;
	      }
	    }//end loop over transformed obs

		weights[i]          = weight;
		particles[i].weight = weight;

	}//End loop over all particles
}

void ParticleFilter::resample()
{
	// Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::vector<Particle> Particles_Resampled;
	//The discrete distribution allows us to create a distribution
	std::discrete_distribution<int> distribution_weights(weights.begin(), weights.end());
	std::default_random_engine gen;

	int total_particles = particles.size();
	for(int i = 0; i < total_particles; i++)
	{
	  Particles_Resampled.push_back(particles[distribution_weights(gen)]);
	}//end loop over total particles

	particles = Particles_Resampled;
}

void ParticleFilter::write(std::string filename)
{
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
