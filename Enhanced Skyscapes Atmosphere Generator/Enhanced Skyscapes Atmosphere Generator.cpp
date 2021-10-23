#include <glm/glm.hpp>

#include <iostream>
#include <fstream>

#include <vector>
#include <execution>

const double earth_radius = 6371008.7714;
const glm::dvec3 earth_center = glm::dvec3(0.0, -1.0 * earth_radius, 0.0);

const double atmosphere_height = 50000.0;

const double atmosphere_base = -0.01;
const double atmosphere_top = atmosphere_height + 0.01;

const glm::dvec3 rayleigh_coefficient = glm::dvec3(5.8e-6, 1.35e-5, 3.31e-5);
const double mie_coefficient = 2.1e-5;
const glm::dvec3 ozone_coefficient = 0.06 * glm::dvec3(3.426e-5, 8.298e-5, 0.356e-5);

const double rayleigh_scale_height = 8000.0;
const double mie_scale_height = 1200.0;

const double ozone_base = 25000.0;
const double ozone_scale_height = 3250.0;

const int step_count = 1000;

double get_ray_height(glm::dvec3 ray_position)
{
	return glm::length(ray_position - earth_center) - earth_radius;
}

double map(double input_value, double input_start, double input_end, double output_start, double output_end)
{
	double slope = (output_end - output_start) / (input_end - input_start);

	return glm::clamp(output_start + (slope * (input_value - input_start)), glm::min(output_start, output_end), glm::max(output_start, output_end));
}

glm::dvec2 ray_sphere_intersections(glm::dvec3 ray_position, glm::dvec3 ray_direction, double sphere_height)
{
	glm::dvec3 ray_earth_vector = ray_position - earth_center;

	double coefficient_1 = 2.0 * glm::dot(ray_direction, ray_earth_vector);
	double coefficient_2 = glm::dot(ray_earth_vector, ray_earth_vector) - glm::pow(earth_radius + sphere_height, 2.0);

	double discriminant = glm::pow(coefficient_1, 2.0) - (4.0 * coefficient_2);

	if (discriminant < 0.0) return glm::dvec2(0.0, 0.0);
	else
	{
		double lower_solution = ((-1.0 * coefficient_1) - glm::sqrt(discriminant)) / 2.0;
		double higher_solution = ((-1.0 * coefficient_1) + glm::sqrt(discriminant)) / 2.0;

		if (lower_solution < 0.0) return glm::dvec2(glm::max(higher_solution, 0.0), 0.0);
		else return glm::dvec2(lower_solution, higher_solution);
	}
}

double ray_atmosphere_intersection(glm::dvec3 ray_position, glm::dvec3 ray_direction)
{
	glm::dvec2 inner_sphere_intersections = ray_sphere_intersections(ray_position, ray_direction, atmosphere_base);
	glm::dvec2 outer_sphere_intersections = ray_sphere_intersections(ray_position, ray_direction, atmosphere_top);

	double lower_distance = glm::min(inner_sphere_intersections.x, outer_sphere_intersections.x);
	double higher_distance = glm::max(inner_sphere_intersections.x, outer_sphere_intersections.x);

	if (lower_distance == 0.0) return higher_distance;
	else return lower_distance;
}

double atmosphere_density(double ray_height, double scale_height)
{
	return exp(-1.0 * (glm::max(ray_height, 0.0) / scale_height));
}

double ozone_density(double ray_height)
{
	return exp(-1.0 * (glm::abs(ray_height - ozone_base) / ozone_scale_height));
}

std::tuple<glm::dvec3, glm::dvec3> render_atmosphere(glm::dvec3 view_position, glm::dvec3 view_direction, glm::dvec3 sun_direction)
{
	glm::dvec3 view_ray_position = view_position;
	double view_step_size = ray_atmosphere_intersection(view_position, view_direction) / double(step_count);

	double view_rayleigh_depth = 0.0;
	double view_mie_depth = 0.0;
	double view_ozone_depth = 0.0;

	glm::dvec3 rayleigh_color = glm::dvec3(0.0, 0.0, 0.0);
	glm::dvec3 mie_color = glm::dvec3(0.0, 0.0, 0.0);

	double previous_view_ray_height = get_ray_height(view_ray_position);

	double previous_view_rayleigh_density = atmosphere_density(previous_view_ray_height, rayleigh_scale_height);
	double previous_view_mie_density = atmosphere_density(previous_view_ray_height, mie_scale_height);
	double previous_view_ozone_density = ozone_density(previous_view_ray_height);

	for (int view_step_index = 0; view_step_index < step_count; view_step_index++)
	{
		double current_view_ray_height = get_ray_height(view_ray_position + (view_direction * view_step_size));

		double current_view_rayleigh_density = atmosphere_density(current_view_ray_height, rayleigh_scale_height);
		double current_view_mie_density = atmosphere_density(current_view_ray_height, mie_scale_height);
		double current_view_ozone_density = ozone_density(current_view_ray_height);

		double rayleigh_depth = 0.5 * (previous_view_rayleigh_density + current_view_rayleigh_density) * view_step_size;
		double mie_depth = 0.5 * (previous_view_mie_density + current_view_mie_density) * view_step_size;
		double ozone_depth = 0.5 * (previous_view_ozone_density + current_view_ozone_density) * view_step_size;

		view_rayleigh_depth += rayleigh_depth;
		view_mie_depth += mie_depth;
		view_ozone_depth += ozone_depth;

		glm::dvec3 sun_ray_position = view_ray_position;

		if (ray_sphere_intersections(sun_ray_position, sun_direction, atmosphere_base).x == 0.0)
		{
			double sun_rayleigh_depth = 0.0;
			double sun_mie_depth = 0.0;
			double sun_ozone_depth = 0.0;

			double sun_step_size = ray_sphere_intersections(sun_ray_position, sun_direction, atmosphere_top).x / double(step_count);

			double previous_sun_ray_height = get_ray_height(sun_ray_position);

			double previous_sun_rayleigh_density = atmosphere_density(previous_sun_ray_height, rayleigh_scale_height);
			double previous_sun_mie_density = atmosphere_density(previous_sun_ray_height, mie_scale_height);
			double previous_sun_ozone_density = ozone_density(previous_sun_ray_height);

			for (int sun_step_index = 0; sun_step_index < step_count; sun_step_index++)
			{
				double current_sun_ray_height = get_ray_height(sun_ray_position + (sun_direction * sun_step_size));

				double current_sun_rayleigh_density = atmosphere_density(current_sun_ray_height, rayleigh_scale_height);
				double current_sun_mie_density = atmosphere_density(current_sun_ray_height, mie_scale_height);
				double current_sun_ozone_density = ozone_density(current_sun_ray_height);

				sun_rayleigh_depth += 0.5 * (previous_sun_rayleigh_density + current_sun_rayleigh_density) * sun_step_size;
				sun_mie_depth += 0.5 * (previous_sun_mie_density + current_sun_mie_density) * sun_step_size;
				sun_ozone_depth += 0.5 * (previous_sun_ozone_density + current_sun_ozone_density) * sun_step_size;

				previous_sun_rayleigh_density = current_sun_rayleigh_density;
				previous_sun_mie_density = current_sun_mie_density;
				previous_sun_ozone_density = current_sun_ozone_density;

				sun_ray_position += sun_direction * sun_step_size;
			}
		
			glm::dvec3 scattering_integrand = (rayleigh_coefficient * (view_rayleigh_depth + sun_rayleigh_depth)) + ((mie_coefficient / 0.9) * (view_mie_depth + sun_mie_depth)) + (ozone_coefficient * (view_ozone_depth + sun_ozone_depth));
			glm::dvec3 scattered_light = glm::exp(-1.0 * scattering_integrand);

			rayleigh_color += scattered_light * rayleigh_depth;
			mie_color += scattered_light * mie_depth;
		}

		previous_view_rayleigh_density = current_view_rayleigh_density;
		previous_view_mie_density = current_view_mie_density;
		previous_view_ozone_density = current_view_ozone_density;

		view_ray_position += view_direction * view_step_size;
	}

	return std::make_tuple(rayleigh_color, mie_color);
}

glm::dvec3 render_transmittance(glm::dvec3 view_position, glm::dvec3 view_direction)
{
	glm::dvec3 view_ray_position = view_position;
	double view_step_size = ray_atmosphere_intersection(view_position, view_direction) / double(step_count);

	double view_rayleigh_depth = 0.0;
	double view_mie_depth = 0.0;
	double view_ozone_depth = 0.0;

	double previous_view_ray_height = get_ray_height(view_ray_position);

	double previous_view_rayleigh_density = atmosphere_density(previous_view_ray_height, rayleigh_scale_height);
	double previous_view_mie_density = atmosphere_density(previous_view_ray_height, mie_scale_height);
	double previous_view_ozone_density = ozone_density(previous_view_ray_height);

	for (int view_step_index = 0; view_step_index < step_count; view_step_index++)
	{
		double current_view_ray_height = get_ray_height(view_ray_position + (view_direction * view_step_size));

		double current_view_rayleigh_density = atmosphere_density(current_view_ray_height, rayleigh_scale_height);
		double current_view_mie_density = atmosphere_density(current_view_ray_height, mie_scale_height);
		double current_view_ozone_density = ozone_density(current_view_ray_height);

		view_rayleigh_depth += 0.5 * (previous_view_rayleigh_density + current_view_rayleigh_density) * view_step_size;
		view_mie_depth += 0.5 * (previous_view_mie_density + current_view_mie_density) * view_step_size;
		view_ozone_depth += 0.5 * (previous_view_ozone_density + current_view_ozone_density) * view_step_size;

		previous_view_rayleigh_density = current_view_rayleigh_density;
		previous_view_mie_density = current_view_mie_density;
		previous_view_ozone_density = current_view_ozone_density;

		view_ray_position += view_direction * view_step_size;
	}

	glm::dvec3 transmittance_integrand = (rayleigh_coefficient * view_rayleigh_depth) + ((mie_coefficient / 0.9) * view_mie_depth) + (ozone_coefficient * view_ozone_depth);

	return glm::exp(-1.0 * transmittance_integrand);
}

float rayleigh_image[128][128][128][3];
float mie_image[128][128][128][3];

float transmittance_image[1024][1024][3];

const unsigned int rayleigh_file_header[4] = {128, 128, 128, 3};
const unsigned int mie_file_header[4] = {128, 128, 128, 3};

const unsigned int transmittance_file_header[3] = {1024, 1024, 3};

int main()
{
    std::cout << "Generating files!" << std::endl;

	std::vector<int> view_height_indices = std::vector<int>(128);
	std::iota(view_height_indices.begin(), view_height_indices.end(), 0);

	std::for_each(std::execution::par_unseq, view_height_indices.begin(), view_height_indices.end(), [&](int view_height_index)
	{
		for (int view_angle_index = 0; view_angle_index < 128; view_angle_index++)
		{
			for (int sun_angle_index = 0; sun_angle_index < 128; sun_angle_index++)
			{
				double view_height_coordinate = double(view_height_index) / 127.0;
				double view_angle_coordinate = double(view_angle_index) / 127.0;
				double sun_angle_coordinate = double(sun_angle_index) / 127.0;

				double view_height = atmosphere_height * glm::pow(view_height_coordinate, 2.0);

				double view_angle = (2.0 * view_angle_coordinate) - 1.0;
				view_angle = glm::acos(glm::sign(view_angle) * glm::pow(view_angle, 2.0));

				double sun_angle = (2.0 * sun_angle_coordinate) - 1.0;
				sun_angle = glm::acos(glm::sign(sun_angle) * glm::pow(sun_angle, 2.0));

				glm::dvec3 view_position = glm::dvec3(0.0, view_height, 0.0);
				glm::dvec3 view_direction = glm::dvec3(0.0, glm::cos(view_angle), glm::sin(view_angle));
				glm::dvec3 sun_direction = glm::dvec3(0.0, glm::cos(sun_angle), glm::sin(sun_angle));

				glm::dvec3 rayleigh_color;
				glm::dvec3 mie_color;

				std::tie(rayleigh_color, mie_color) = render_atmosphere(view_position, view_direction, sun_direction);

				rayleigh_image[sun_angle_index][view_angle_index][view_height_index][0] = rayleigh_color.x;
				rayleigh_image[sun_angle_index][view_angle_index][view_height_index][1] = rayleigh_color.y;
				rayleigh_image[sun_angle_index][view_angle_index][view_height_index][2] = rayleigh_color.z;

				mie_image[sun_angle_index][view_angle_index][view_height_index][0] = mie_color.x;
				mie_image[sun_angle_index][view_angle_index][view_height_index][1] = mie_color.y;
				mie_image[sun_angle_index][view_angle_index][view_height_index][2] = mie_color.z;
			}
		}
	});

	view_height_indices = std::vector<int>(1024);
	std::iota(view_height_indices.begin(), view_height_indices.end(), 0);

	std::for_each(std::execution::par_unseq, view_height_indices.begin(), view_height_indices.end(), [&](int view_height_index)
	{
		for (int view_angle_index = 0; view_angle_index < 1024; view_angle_index++)
		{
			double view_height_coordinate = double(view_height_index) / 1023.0;
			double view_angle_coordinate = double(view_angle_index) / 1023.0;

			double view_height = atmosphere_height * glm::pow(view_height_coordinate, 2.0);

			double view_angle = (2.0 * view_angle_coordinate) - 1.0;
			view_angle = glm::acos(glm::sign(view_angle) * glm::pow(view_angle, 2.0));

			glm::dvec3 view_position = glm::dvec3(0.0, view_height, 0.0);
			glm::dvec3 view_direction = glm::dvec3(0.0, glm::cos(view_angle), glm::sin(view_angle));

			glm::dvec3 transmittance = render_transmittance(view_position, view_direction);

			transmittance_image[view_angle_index][view_height_index][0] = transmittance.x;
			transmittance_image[view_angle_index][view_height_index][1] = transmittance.y;
			transmittance_image[view_angle_index][view_height_index][2] = transmittance.z;
		}
	});

	std::ofstream rayleigh_file;
	std::ofstream mie_file;

	std::ofstream transmittance_file;

	rayleigh_file.open("rayleigh.bin", std::ofstream::binary);
	mie_file.open("mie.bin", std::ofstream::binary);

	transmittance_file.open("transmittance.bin", std::ofstream::binary);

	rayleigh_file.write((const char*)(rayleigh_file_header), sizeof(rayleigh_file_header));
	rayleigh_file.write((const char*)(rayleigh_image), sizeof(rayleigh_image));

	mie_file.write((const char*)(mie_file_header), sizeof(mie_file_header));
	mie_file.write((const char*)(mie_image), sizeof(mie_image));

	transmittance_file.write((const char*)(transmittance_file_header), sizeof(transmittance_file_header));
	transmittance_file.write((const char*)(transmittance_image), sizeof(transmittance_image));

	rayleigh_file.close();
	mie_file.close();

	transmittance_file.close();

	std::cout << "Done!" << std::endl;

	return 0;
}