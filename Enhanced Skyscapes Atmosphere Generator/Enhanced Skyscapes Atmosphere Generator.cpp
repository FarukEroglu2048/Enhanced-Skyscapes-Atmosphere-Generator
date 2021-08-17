#include <glm/glm.hpp>

#include <iostream>
#include <fstream>

#include <vector>
#include <execution>

const float pi = 3.141592653589793f;

const float earth_radius = 6378100.0f;
const glm::vec3 earth_center = glm::vec3(0.0f, -1.0f * earth_radius, 0.0f);

const float atmosphere_height = 50000.0f;

const glm::vec3 rayleigh_coefficient = glm::vec3(5.8e-6f, 1.35e-5f, 3.31e-5f);
const float mie_coefficient = 2.1e-5f;
const glm::vec3 ozone_coefficient = 0.06f * glm::vec3(3.426e-5f, 8.298e-5f, 0.356e-5f);

const float rayleigh_scale_height = 8000.0f;
const float mie_scale_height = 1200.0f;

const float ozone_base = 25000.0f;
const float ozone_scale_height = 7500.0f;

const int step_count = 50;

float map(float input_value, float input_start, float input_end, float output_start, float output_end)
{
	float slope = (output_end - output_start) / (input_end - input_start);

	return glm::clamp(output_start + (slope * (input_value - input_start)), glm::min(output_start, output_end), glm::max(output_start, output_end));
}

float get_ray_height(glm::vec3 ray_position)
{
	return glm::length(ray_position - earth_center) - earth_radius;
}

glm::vec2 ray_sphere_intersections(glm::vec3 ray_position, glm::vec3 ray_direction, float sphere_height)
{
	glm::vec3 ray_earth_vector = ray_position - earth_center;

	float coefficient_1 = 2.0f * glm::dot(ray_direction, ray_earth_vector);
	float coefficient_2 = glm::dot(ray_earth_vector, ray_earth_vector) - glm::pow(earth_radius + sphere_height, 2.0f);

	float discriminant = glm::pow(coefficient_1, 2.0f) - (4.0f * coefficient_2);

	if (discriminant < 0.0f) return glm::vec2(0.0f, 0.0f);
	else
	{
		float lower_solution = ((-1.0f * coefficient_1) - glm::sqrt(discriminant)) / 2.0f;
		float higher_solution = ((-1.0f * coefficient_1) + glm::sqrt(discriminant)) / 2.0f;

		if (lower_solution < 0.0f) return glm::vec2(glm::max(higher_solution, 0.0f), 0.0f);
		else return glm::vec2(lower_solution, higher_solution);
	}
}

glm::vec2 ray_atmosphere_intersections(glm::vec3 ray_position, glm::vec3 ray_direction)
{
	glm::vec2 layer_intersections = glm::vec2(0.0f, 0.0f);

	glm::vec2 inner_sphere_intersections = ray_sphere_intersections(ray_position, ray_direction, 0.0f);
	glm::vec2 outer_sphere_intersections = ray_sphere_intersections(ray_position, ray_direction, atmosphere_height);

	float ray_height = get_ray_height(ray_position);

	if (ray_height < 0.0f)
	{
		layer_intersections.x = inner_sphere_intersections.x;
		layer_intersections.y = outer_sphere_intersections.x - inner_sphere_intersections.x;
	}
	else if ((ray_height >= 0.0f) && (ray_height <= atmosphere_height))
	{
		float lower_distance = glm::min(inner_sphere_intersections.x, outer_sphere_intersections.x);
		float higher_distance = glm::max(inner_sphere_intersections.x, outer_sphere_intersections.x);

		if (lower_distance == 0.0f) layer_intersections.y = higher_distance;
		else layer_intersections.y = lower_distance;
	}
	else if (ray_height > atmosphere_height)
	{
		layer_intersections.x = outer_sphere_intersections.x;

		if (inner_sphere_intersections.x == 0.0f) layer_intersections.y = outer_sphere_intersections.y - outer_sphere_intersections.x;
		else layer_intersections.y = inner_sphere_intersections.x - outer_sphere_intersections.x;
	}

	return layer_intersections;
}

float atmosphere_density(float ray_height, float scale_height)
{
	return exp(-1.0f * (glm::max(ray_height, 0.0f) / scale_height));
}

float ozone_density(float ray_height)
{
	return exp(-1.0f * (glm::abs(ray_height - ozone_base) / ozone_scale_height));
}

std::tuple<glm::vec3, glm::vec3> render_atmosphere(glm::vec3 view_position, glm::vec3 view_direction, glm::vec3 sun_direction)
{
	glm::vec2 view_ray_intersections = ray_atmosphere_intersections(view_position, view_direction);

	glm::vec3 view_ray_position = view_position + (view_direction * view_ray_intersections.x);
	float view_step_size = view_ray_intersections.y / float(step_count);

	float view_rayleigh_depth = 0.0f;
	float view_mie_depth = 0.0f;
	float view_ozone_depth = 0.0f;

	glm::vec3 rayleigh_color = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 mie_color = glm::vec3(0.0f, 0.0f, 0.0f);

	float previous_view_ray_height = get_ray_height(view_ray_position);

	float previous_view_rayleigh_density = atmosphere_density(previous_view_ray_height, rayleigh_scale_height);
	float previous_view_mie_density = atmosphere_density(previous_view_ray_height, mie_scale_height);
	float previous_view_ozone_density = ozone_density(previous_view_ray_height);

	for (int view_step_index = 0; view_step_index < step_count; view_step_index++)
	{
		float current_view_ray_height = get_ray_height(view_ray_position + (view_direction * view_step_size));

		float current_view_rayleigh_density = atmosphere_density(current_view_ray_height, rayleigh_scale_height);
		float current_view_mie_density = atmosphere_density(current_view_ray_height, mie_scale_height);
		float current_view_ozone_density = ozone_density(current_view_ray_height);

		float rayleigh_depth = 0.5f * (previous_view_rayleigh_density + current_view_rayleigh_density) * view_step_size;
		float mie_depth = 0.5f * (previous_view_mie_density + current_view_mie_density) * view_step_size;
		float ozone_depth = 0.5f * (previous_view_ozone_density + current_view_ozone_density) * view_step_size;

		view_rayleigh_depth += rayleigh_depth;
		view_mie_depth += mie_depth;
		view_ozone_depth += ozone_depth;

		glm::vec3 sun_ray_position = view_position;
		float sun_step_size = ray_sphere_intersections(sun_ray_position, sun_direction, atmosphere_height).x / float(step_count);

		float sun_rayleigh_depth = 0.0f;
		float sun_mie_depth = 0.0f;
		float sun_ozone_depth = 0.0f;

		float previous_sun_ray_height = get_ray_height(sun_ray_position);

		float previous_sun_rayleigh_density = atmosphere_density(previous_sun_ray_height, rayleigh_scale_height);
		float previous_sun_mie_density = atmosphere_density(previous_sun_ray_height, mie_scale_height);
		float previous_sun_ozone_density = ozone_density(previous_sun_ray_height);

		for (int sun_step_index = 0; sun_step_index < step_count; sun_step_index++)
		{
			float current_sun_ray_height = get_ray_height(sun_ray_position + (sun_direction * sun_step_size));

			float current_sun_rayleigh_density = atmosphere_density(current_sun_ray_height, rayleigh_scale_height);
			float current_sun_mie_density = atmosphere_density(current_sun_ray_height, mie_scale_height);
			float current_sun_ozone_density = ozone_density(current_sun_ray_height);

			sun_rayleigh_depth += 0.5f * (previous_sun_rayleigh_density + current_sun_rayleigh_density) * sun_step_size;
			sun_mie_depth += 0.5f * (previous_sun_mie_density + current_sun_mie_density) * sun_step_size;
			sun_ozone_depth += 0.5f * (previous_sun_ozone_density + current_sun_ozone_density) * sun_step_size;

			previous_sun_rayleigh_density = current_sun_rayleigh_density;
			previous_sun_mie_density = current_sun_mie_density;
			previous_sun_ozone_density = current_sun_ozone_density;

			sun_ray_position += sun_direction * sun_step_size;
		}

		glm::vec3 scattering_integrand = (rayleigh_coefficient * (view_rayleigh_depth + sun_rayleigh_depth)) + ((mie_coefficient / 0.9f) * (view_mie_depth + sun_mie_depth)) + (ozone_coefficient * (view_ozone_depth + sun_ozone_depth));
		glm::vec3 scattered_light = glm::exp(-1.0f * scattering_integrand);

		rayleigh_color += scattered_light * rayleigh_depth;
		mie_color += scattered_light * mie_depth;

		previous_view_rayleigh_density = current_view_rayleigh_density;
		previous_view_mie_density = current_view_mie_density;
		previous_view_ozone_density = current_view_ozone_density;

		view_ray_position += view_direction * view_step_size;
	}

	return std::make_tuple(rayleigh_color, mie_color);
}

glm::vec3 render_transmittance(glm::vec3 view_position, glm::vec3 view_direction)
{
	glm::vec2 view_ray_intersections = ray_atmosphere_intersections(view_position, view_direction);

	glm::vec3 view_ray_position = view_position + (view_direction * view_ray_intersections.x);
	float view_step_size = view_ray_intersections.y / float(step_count);

	float view_rayleigh_depth = 0.0f;
	float view_mie_depth = 0.0f;
	float view_ozone_depth = 0.0f;

	glm::vec3 rayleigh_color = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 mie_color = glm::vec3(0.0f, 0.0f, 0.0f);

	float previous_view_ray_height = get_ray_height(view_ray_position);

	float previous_view_rayleigh_density = atmosphere_density(previous_view_ray_height, rayleigh_scale_height);
	float previous_view_mie_density = atmosphere_density(previous_view_ray_height, mie_scale_height);
	float previous_view_ozone_density = ozone_density(previous_view_ray_height);

	for (int view_step_index = 0; view_step_index < step_count; view_step_index++)
	{
		float current_view_ray_height = get_ray_height(view_ray_position + (view_direction * view_step_size));

		float current_view_rayleigh_density = atmosphere_density(current_view_ray_height, rayleigh_scale_height);
		float current_view_mie_density = atmosphere_density(current_view_ray_height, mie_scale_height);
		float current_view_ozone_density = ozone_density(current_view_ray_height);

		view_rayleigh_depth += 0.5f * (previous_view_rayleigh_density + current_view_rayleigh_density) * view_step_size;
		view_mie_depth += 0.5f * (previous_view_mie_density + current_view_mie_density) * view_step_size;
		view_ozone_depth += 0.5f * (previous_view_ozone_density + current_view_ozone_density) * view_step_size;

		previous_view_rayleigh_density = current_view_rayleigh_density;
		previous_view_mie_density = current_view_mie_density;
		previous_view_ozone_density = current_view_ozone_density;

		view_ray_position += view_direction * view_step_size;
	}

	glm::vec3 transmittance_integrand = (rayleigh_coefficient * view_rayleigh_depth) + ((mie_coefficient / 0.9f) * view_mie_depth) + (ozone_coefficient * view_ozone_depth);

	return glm::exp(-1.0f * transmittance_integrand);
}

float rayleigh_image[128][128][128][3];
float mie_image[128][128][128][3];

unsigned char transmittance_image[1024][1024][3];

const unsigned int rayleigh_file_header[3] = {128, 128, 128};
const unsigned int mie_file_header[3] = {128, 128, 128};

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
				float view_height_coordinate = float(view_height_index) / 127.0;
				float view_angle_coordinate = float(view_angle_index) / 127.0;
				float sun_angle_coordinate = float(sun_angle_index) / 127.0;

				float view_height = atmosphere_height * glm::pow(view_height_coordinate, 2.0);

				float view_angle = (2.0 * view_angle_coordinate) - 1.0;
				view_angle = 0.5 * pi * glm::sign(view_angle) * glm::pow(view_angle, 2.0);

				float sun_angle = (2.0 * sun_angle_coordinate) - 1.0;
				sun_angle = 0.5 * pi * glm::sign(sun_angle) * glm::pow(sun_angle, 2.0);

				glm::vec3 view_position = glm::vec3(0.0, glm::max(view_height, 1.0f), 0.0);
				glm::vec3 view_direction = glm::vec3(0.0, glm::sin(view_angle), glm::cos(view_angle));
				glm::vec3 sun_direction = glm::vec3(0.0, glm::sin(sun_angle), glm::cos(sun_angle));

				glm::vec3 rayleigh_color;
				glm::vec3 mie_color;

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
			float view_height_coordinate = float(view_height_index) / 1023.0;
			float view_angle_coordinate = float(view_angle_index) / 1023.0;

			float view_height = atmosphere_height * view_height_coordinate;
			float view_angle = 0.5 * pi * ((2.0 * view_angle_coordinate) - 1.0);

			glm::vec3 view_position = glm::vec3(0.0, glm::max(view_height, 1.0f), 0.0);
			glm::vec3 view_direction = glm::vec3(0.0, glm::sin(view_angle), glm::cos(view_angle));

			glm::vec3 transmittance = render_transmittance(view_position, view_direction);

			transmittance_image[view_angle_index][view_height_index][0] = (unsigned char)(255.0 * transmittance.x);
			transmittance_image[view_angle_index][view_height_index][1] = (unsigned char)(255.0 * transmittance.y);
			transmittance_image[view_angle_index][view_height_index][2] = (unsigned char)(255.0 * transmittance.z);
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

	transmittance_file.write((const char*)(transmittance_image), sizeof(transmittance_image));

	rayleigh_file.close();
	mie_file.close();

	transmittance_file.close();

	return 0;
}