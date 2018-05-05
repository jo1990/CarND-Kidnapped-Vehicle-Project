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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles=100;
	default_random_engine gen;

	normal_distribution<double> N_x(x,std[0]);
	normal_distribution<double> N_y(y,std[1]);
	normal_distribution<double> N_theta(theta,std[2]);



	for(int i=0;i<num_particles;i++)
	{
		Particle particle;

		particle.id=i;
//		particle.x=x;
//		particle.y=y;
//		particle.theta=theta;
//		particle.weight=1.0;


		particle.x+=N_x(gen);
		particle.y+=N_y(gen);
		particle.theta+=N_theta(gen);

		particles.push_back(particle);
		//		weights.push_back(1.0);

	}

	is_initialized=true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;


	for(int i=0;i<num_particles;i++)
	{
		//		double new_x;
		//		double new_y;
		//		double new_theta;


		normal_distribution<double> new_x(0, std_pos[0]);
		normal_distribution<double> new_y(0, std_pos[1]);
		normal_distribution<double> new_theta(0, std_pos[2]);


		// calculate new state
		if (yaw_rate== 0) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		else {
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}

		// add noise
		particles[i].x += new_x(gen);
		particles[i].y += new_y(gen);
		particles[i].theta += new_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	int num_predictions=predicted.size();
	int num_obsevations=observations.size();

	for(int i=0;i<num_obsevations;i++)
	{
		int landMark_id=-1;
		double minDist_LmPred_LmObs=numeric_limits<double>::max();

		LandmarkObs LandMark_obs=observations[i];

		for(int p=0;p<num_predictions;p++)
		{
			LandmarkObs LandMark_pred=predicted[p];

			double meas_dist=dist(LandMark_obs.x,LandMark_obs.y,LandMark_pred.x,LandMark_pred.y);

			if(meas_dist<minDist_LmPred_LmObs)
			{
				minDist_LmPred_LmObs=meas_dist;
				landMark_id=LandMark_pred.id;
			}
		}

		observations[i].id=landMark_id;

	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html





	for(int i=0;i<num_particles;i++)
	{
		vector<LandmarkObs> landMark_pred;

		for(int j=0;j<map_landmarks.landmark_list.size();j++)
		{
			double LandMark_x=map_landmarks.landmark_list[j].x_f;
			double LandMark_y=map_landmarks.landmark_list[j].y_f;
			int LandMark_id=map_landmarks.landmark_list[j].id_i;

			double meas_dist=dist(LandMark_x,LandMark_y,particles[i].x,particles[i].y);

			if(meas_dist <= sensor_range)
			{
				landMark_pred.push_back(LandmarkObs{LandMark_id,LandMark_x,LandMark_y});
			}


		}

		vector<LandmarkObs> trans_landMarks;

		//xm=xp​+(cosθ × Xc)−(sinθ × Yc)
		//ym=yp​+(sinθ × Xc)+(cosθ × Yc)
		for(int t=0;t<observations.size();t++)
		{
			double trans_x=particles[i].x+(cos(particles[i].theta ) * observations[t].x - sin(particles[i].theta) * observations[t].y);

			double trans_y=particles[i].y+(sin(particles[i].theta ) * observations[t].x + cos(particles[i].theta) * observations[t].y);

			trans_landMarks.push_back(LandmarkObs{observations[t].id,trans_x,trans_y});
		}

		dataAssociation(landMark_pred,trans_landMarks);

		particles[i].weight=1.0;

		for(int l=0;l<trans_landMarks.size();l++)
		{
			LandmarkObs lm_pred;

			int pred_id=trans_landMarks[l].id;

			for(int k=0;k<landMark_pred.size();k++)
			{
				if(landMark_pred[k].id==pred_id)
				{
					lm_pred.x=landMark_pred[k].x;
					lm_pred.y=landMark_pred[k].y;

				}
			}


			//			# calculate normalization term
			//			gauss_norm= (1/(2 * np.pi * sig_x * sig_y))
			//
			//			# calculate exponent
			//			exponent= ((x_obs - mu_x)**2)/(2 * sig_x**2) + ((y_obs - mu_y)**2)/(2 * sig_y**2)
			//
			//			# calculate weight using normalization terms and exponent
			//			weight= gauss_norm * math.exp(-exponent)


			double gauss_norm=1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);

			double exponent=exp( -( pow( lm_pred.x - trans_landMarks[l].x, 2 ) / ( 2 * pow( std_landmark[0], 2 ) ) +
					( pow( lm_pred.y - trans_landMarks[l].y, 2 ) / ( 2  *pow( std_landmark[1], 2 ) ) ) ) );

			double weight=gauss_norm * exponent;

			particles[i].weight  = particles[i].weight * weight;



		}
	}


}

//modified from the interent
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;

	vector<Particle> new_particles;

	// get all of the current weights
	vector<double> weights;
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}

	// generate random starting index for resampling wheel
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	auto index = uniintdist(gen);

	// get max weight
	double max_weight = *max_element(weights.begin(), weights.end());

	// uniform random distribution [0.0, max_weight)
	uniform_real_distribution<double> unirealdist(0.0, max_weight);

	double beta = 0.0;

	// spin the resample wheel!
	for (int i = 0; i < num_particles; i++) {
		beta += unirealdist(gen) * 2.0;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		new_particles.push_back(particles[index]);
	}

	particles = new_particles;
}



Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
		const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations= associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
