#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/bounding_box_diagonal.h>
#include <igl/slice.h>
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

/* ------- Required For My Solution ------- */
// Structures
typedef struct Resolution { int x, y, z; } Resolution;
typedef struct BoundingBox
{
	float diagonal;
	float diagonal_enlargement;
	Eigen::RowVector3d min;
	Eigen::RowVector3d max;
	Eigen::RowVector3d dimensions;
} BoundingBox;
typedef struct SpatialIndex
{
	std::vector<std::vector<int>> structure;
	Resolution resolution;
	float enlarge_by;
} SpatialIndex;

// Global Variables
std::string model_name = "cat";
BoundingBox bb;
SpatialIndex spatial_index;
float pos_epsilon_coeff = 0.03f;
float neg_epsilon_coeff = 0.01f;

// Functions
void extruct_model_name(string file_name);
void initialize_parameters();
// Queries
std::pair<Eigen::VectorXd, float> find_closest_point(unsigned int i, float epsilon, short sign);
std::pair<std::vector<int>, std::vector<double>> find_points_within_distance(Eigen::RowVector3d q, float h);
/* ------- Required For My Solution END ------- */

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
Resolution resolution = { 20, 20, 20 };

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);



// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
	grid_points.resize(0, 3);
	grid_colors.resize(0, 3);
	grid_lines.resize(0, 6);
	grid_values.resize(0);
	V.resize(0, 3);
	F.resize(0, 3);
	FN.resize(0, 3);

	// Grid spacing
	const double dx = bb.dimensions(0) / (double)(resolution.x - 1);
	const double dy = bb.dimensions(1) / (double)(resolution.y - 1);
	const double dz = bb.dimensions(2) / (double)(resolution.z - 1);
	// 3D positions of the grid points -- see slides or marching_cubes.h for ordering
	grid_points.resize(resolution.x * resolution.y * resolution.z, 3);
	// Create each gridpoint
	for (unsigned int x = 0; x < resolution.x; ++x) {
		for (unsigned int y = 0; y < resolution.y; ++y) {
			for (unsigned int z = 0; z < resolution.z; ++z) {
				// Linear index of the point at (x,y,z)
				int index = x + resolution.y * (y + resolution.z * z);
				// 3D point at (x,y,z)
				grid_points.row(index) = bb.min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
			}
		}
	}
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc() {
	int coeffs = 0;
	switch (polyDegree) {
	case 0: coeffs = 1; break;
	case 1: coeffs = 4; break;
	case 2: coeffs = 10; break;
	default:
		cerr << "Unsupported polynomial degree" << endl;
		return;
	}

	double radius = wendlandRadius * bb.diagonal;

	// Scalar values of the grid points (the implicit function values)
	grid_values.resize(resolution.x * resolution.y * resolution.z);

	// Evaluate signed distance function at each gridpoint
	for (unsigned int x = 0; x < resolution.x; ++x) {
		for (unsigned int y = 0; y < resolution.y; ++y) {
			for (unsigned int z = 0; z < resolution.z; ++z) {
				// Linear index of the point at (x,y,z)
				int index = x + resolution.x * (y + resolution.y * z);
				auto query_result = find_points_within_distance(grid_points.row(index), radius);
				if (query_result.first.size() > 0)
				{
					Eigen::MatrixXd points;
					igl::slice(constrained_points, Eigen::VectorXi::Map(query_result.first.data(), query_result.first.size()), Eigen::Matrix<int, 1, 3>(0, 1, 2), points);
					Eigen::MatrixXd points_reorder(points.rows(), 3);
					points_reorder << points.col(1), points.col(2), points.col(0);
					Eigen::MatrixXd points_square = points.cwiseProduct(points), points_product = points.cwiseProduct(points_reorder);

					Eigen::MatrixXd a(points.rows(), coeffs);
					Eigen::VectorXd b(coeffs);
					double x = grid_points(index, 0), y = grid_points(index, 1), z = grid_points(index, 2);
					switch (polyDegree) {
					case 0:
						a << Eigen::MatrixXd::Ones(points.rows(), 1);
						b << 1;
						break;
					case 1:
						a << Eigen::MatrixXd::Ones(points.rows(), 1), points;
						b << 1, x, y, z;
						break;
					case 2:
						a << Eigen::MatrixXd::Ones(points.rows(), 1), points, points_square, points_product;
						b << 1, x, y, z, pow(x, 2), pow(y, 2), pow(z, 2), (x * y), (y * z), (z * x);
						break;
					}

					auto ratios = (Eigen::VectorXd::Map(query_result.second.data(), query_result.second.size())).array() / radius;
					Eigen::VectorXd wendland_functions = (1 - ratios).pow(4) * (4 * ratios + 1);

					Eigen::VectorXd f;
					igl::slice(constrained_values, Eigen::VectorXi::Map(query_result.first.data(), query_result.first.size()), Eigen::Matrix<int, 1, 1>(0), f);

					Eigen::VectorXd solution = (a.transpose() * wendland_functions.asDiagonal() * a).ldlt().solve(a.transpose() * wendland_functions.asDiagonal() * f);

					// Value at (x,y,z) = implicit function
					grid_values[index] = b.dot(solution);
				}
				else
				{
					// no constraint points within radius
					grid_values[index] = 100;
				}
			}
		}
	}
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
	int nnodes = grid_points.rows();
	grid_lines.resize(3 * nnodes, 6);
	int numLines = 0;

	for (unsigned int x = 0; x < resolution.x; ++x) {
		for (unsigned int y = 0; y < resolution.y; ++y) {
			for (unsigned int z = 0; z < resolution.z; ++z) {
				int index = x + resolution.x * (y + resolution.y * z);
				if (x < resolution.x - 1) {
					int index1 = (x + 1) + y * resolution.x + z * resolution.y * resolution.z;
					grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
				}
				if (y < resolution.y - 1) {
					int index1 = x + (y + 1) * resolution.x + z * resolution.y * resolution.z;
					grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
				}
				if (z < resolution.z - 1) {
					int index1 = x + y * resolution.x + (z + 1) * resolution.y * resolution.z;
					grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
				}
			}
		}
	}

	grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
	if (key == '1') {
		// Show imported points
		viewer.data().clear();
		viewer.core.align_camera_center(P);
		viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
	}

	if (key == '2')
	{
		// Show all constraints
		viewer.data().clear();
		viewer.core.align_camera_center(P);
		// Add your code for computing auxiliary constraint points here
		// Add code for displaying all points, as above

		// implementing a spatial index
		Eigen::RowVector3d enlarged_dim = bb.dimensions / bb.diagonal_enlargement;
		spatial_index.resolution = { (int)enlarged_dim(0) + 1, (int)enlarged_dim(1) + 1, (int)enlarged_dim(2) + 1 };
		spatial_index.structure.clear();
		spatial_index.structure.resize(spatial_index.resolution.x * spatial_index.resolution.y * spatial_index.resolution.z);
		Eigen::MatrixXd spatial = (P - bb.min.replicate(P.rows(), 1)) / bb.diagonal_enlargement;
		for (int i = 0; i < spatial.rows(); i++)
		{
			int index = (int)spatial(i, 0) + spatial_index.resolution.x * int(spatial(i, 1) + spatial_index.resolution.y * (int)spatial(i, 2));
			spatial_index.structure[index].push_back(i);
		}

		// calculating and displaying the constraints
		constrained_points.resize(3 * P.rows(), 3);
		constrained_values.setZero(3 * P.rows());
		for (int i = 0; i < P.rows(); i++)
		{
			// positive constraint - finding the closest outside offset point
			auto query_result = find_closest_point(i, pos_epsilon_coeff * bb.diagonal, 1);
			constrained_points.row(P.rows() + i) = query_result.first;
			constrained_values(P.rows() + i) = query_result.second;
			// negative constraint - finding the closest inside offset point
			query_result = find_closest_point(i, neg_epsilon_coeff * bb.diagonal, -1);
			constrained_points.row(2 * P.rows() + i) = query_result.first;
			constrained_values(2 * P.rows() + i) = -query_result.second;
		}

		viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 1));
		viewer.data().add_points(constrained_points.block(P.rows(), 0, P.rows(), 3), Eigen::RowVector3d(1, 0, 0));
		viewer.data().add_points(constrained_points.block(2 * P.rows(), 0, P.rows(), 3), Eigen::RowVector3d(0, 1, 0));
	}

	if (key == '3') {
		// Show grid points with colored nodes and connected with lines
		viewer.data().clear();
		viewer.core.align_camera_center(P);
		// Add code for creating a grid
		// Add your code for evaluating the implicit function at the grid points
		// Add code for displaying points and lines
		// You can use the following example:

		/*** begin: sphere example, replace (at least partially) with your code ***/
		// Make grid
		createGrid();

		// Evaluate implicit function
		evaluateImplicitFunc();

		// get grid lines
		getLines();

		// Code for coloring and displaying the grid points and lines
		// Assumes that grid_values and grid_points have been correctly assigned.
		grid_colors.setZero(grid_points.rows(), 3);

		// Build color map
		for (int i = 0; i < grid_points.rows(); ++i) {
			double value = grid_values(i);
			if (value < 0) {
				grid_colors(i, 1) = 1;
			}
			else {
				if (value > 0)
					grid_colors(i, 0) = 1;
			}
		}

		// Draw lines and points
		viewer.data().add_points(grid_points, grid_colors);
		viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
			grid_lines.block(0, 3, grid_lines.rows(), 3),
			Eigen::RowVector3d(0.8, 0.8, 0.8));
		/*** end: sphere example ***/
	}

	if (key == '4') {
		// Show reconstructed mesh
		viewer.data().clear();
		// Code for computing the mesh (V,F) from grid_points and grid_values
		if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
			cerr << "Not enough data for Marching Cubes !" << endl;
			return true;
		}
		// Run marching cubes
		igl::copyleft::marching_cubes(grid_values, grid_points, resolution.x, resolution.y, resolution.z, V, F);
		if (V.rows() == 0) {
			cerr << "Marching Cubes failed!" << endl;
			return true;
		}

		igl::per_face_normals(V, F, FN);
		viewer.data().set_mesh(V, F);
		viewer.data().show_lines = true;
		viewer.data().show_faces = true;
		viewer.data().set_normals(FN);

		igl::writeOFF("../results/" + model_name + ".off", V, F);
	}

	return true;
}

bool callback_load_mesh(Viewer& viewer, string filename)
{
	igl::readOFF(filename, P, F, N);
	extruct_model_name(filename);
	initialize_parameters();
	callback_key_down(viewer, '1', 0);
	return true;
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "Usage ex2_bin <mesh.off>" << endl;
		igl::readOFF("../data/" + model_name + ".off", P, F, N);
	}
	else
	{
		// Read points and normals
		igl::readOFF(argv[1], P, F, N);
		extruct_model_name(argv[1]);
	}
	initialize_parameters();

	Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);
	viewer.data().point_size = 5;
	viewer.callback_key_down = callback_key_down;

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...

			// TODO: Add more parameters to tweak here...
			if (ImGui::TreeNode("Resolution"))
			{
				ImGui::InputInt("x", &resolution.x, 0, 0);
				ImGui::InputInt("y", &resolution.y, 0, 0);
				ImGui::InputInt("z", &resolution.z, 0, 0);
				ImGui::TreePop();
			}
			ImGui::InputDouble("Wendland Radius", &wendlandRadius, 0, 0, "%.3f");
			ImGui::InputInt("Polynomial  Degree", &polyDegree, 0, 0);

			if (ImGui::TreeNode("Epsilon Coefficients"))
			{
				ImGui::InputFloat("Positive", &pos_epsilon_coeff, 0.01f, 0.01f);
				ImGui::InputFloat("Negative", &neg_epsilon_coeff, 0.01f, 0.01f);
				ImGui::TreePop();
			}
			ImGui::InputFloat("Spatial Enlargement", &spatial_index.enlarge_by, 0.1f, 0.1f);

			float w = ImGui::GetContentRegionAvailWidth(), p = ImGui::GetStyle().FramePadding.x;
			if (ImGui::Button("Reset Grid", ImVec2((w - p) / 2.f, 0)))
			{
				std::cout << "ResetGrid\n";
				// Recreate the grid
				createGrid();
				// Switch view to show the grid
				callback_key_down(viewer, '3', 0);
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Load Cloud", ImVec2((w - p) / 2.f, 0)))
			{
				auto file_name = igl::file_dialog_open();
				if (file_name.length() > 0)
				{
					callback_load_mesh(viewer, file_name);
				}
			}
		}

	};

	callback_key_down(viewer, '1', 0);

	viewer.launch();
}

void extruct_model_name(string file_name)
{
	const auto index = file_name.find_last_of('.');
	if (index == std::string::npos)
	{
		model_name = file_name;
	}
	else
	{
		model_name = file_name.substr(0, index);
	}
}

void initialize_parameters()
{
	spatial_index.structure.clear();
	spatial_index.resolution = { 0, 0, 0 };
	spatial_index.enlarge_by = 0.1f;
	bb.max = P.colwise().maxCoeff();
	bb.min = P.colwise().minCoeff();
	bb.dimensions = bb.max - bb.min;
	bb.diagonal = igl::bounding_box_diagonal(P);
	bb.diagonal_enlargement = spatial_index.enlarge_by * bb.diagonal;
}

std::pair<Eigen::VectorXd, float> find_closest_point(unsigned int i, float epsilon, short sign)
{
	bool found = false;
	Eigen::RowVectorXd closest = P.row(i) + sign * epsilon * N.row(i);
	while (!found)
	{
		int closest_i = -1;
		float min_distance = bb.diagonal;
		Eigen::RowVector3d dim = (closest - bb.min) / bb.diagonal_enlargement;
		for (int i = max(0, (int)dim(0) - 1); i < min(spatial_index.resolution.x, (int)dim(0) + 2); i++)
		{
			for (int j = max(0, (int)dim(1) - 1); j < min(spatial_index.resolution.y, (int)dim(1) + 2); j++)
			{
				for (int k = max(0, (int)dim(2) - 1); k < min(spatial_index.resolution.z, (int)dim(2) + 2); k++)
				{
					auto& s_i = spatial_index.structure[i + spatial_index.resolution.x * (j + spatial_index.resolution.y * k)];
					for (int index = 0; index < s_i.size(); index++)
					{
						auto distance = (P.row(s_i[index]) - closest).norm();
						if (distance < min_distance)
						{
							min_distance = distance;
							closest_i = s_i[index];
						}
					}
				}
			}
		}
		if (closest_i != i)
		{
			epsilon /= 2;
			closest = P.row(i) + sign * epsilon * N.row(i);
		}
		else
		{
			found = true;
		}
	}
	return std::make_pair(closest, epsilon);
}

std::pair<std::vector<int>, std::vector<double>> find_points_within_distance(Eigen::RowVector3d q, float h)
{
	std::vector<int> points;
	std::vector<double> distances;
	Eigen::RowVector3d dim = (q - bb.min) / bb.diagonal_enlargement;
	int cells = ceil(h / bb.diagonal_enlargement);

	for (int i = max(0, (int)dim(0) - cells); i < min(spatial_index.resolution.x, (int)dim(0) + cells + 1); i++)
	{
		for (int j = max(0, (int)dim(1) - cells); j < min(spatial_index.resolution.y, (int)dim(1) + cells + 1); j++)
		{
			for (int k = max(0, (int)dim(2) - cells); k < min(spatial_index.resolution.z, (int)dim(2) + cells + 1); k++)
			{
				double distance;
				auto& s_i = spatial_index.structure[i + spatial_index.resolution.x * (j + spatial_index.resolution.y * k)];
				for (int index = 0; index < s_i.size(); index++)
				{
					// zero
					distance = (constrained_points.row(s_i[index]) - q).norm();
					if (distance <= h)
					{
						points.push_back(s_i[index]);
						distances.push_back(distance);
					}
					// positive
					distance = (constrained_points.row(P.rows() + s_i[index]) - q).norm();
					if (distance <= h)
					{
						points.push_back(P.rows() + s_i[index]);
						distances.push_back(distance);
					}
					// negative
					distance = (constrained_points.row(P.rows() * 2 + s_i[index]) - q).norm();
					if (distance <= h)
					{
						points.push_back(P.rows() * 2 + s_i[index]);
						distances.push_back(distance);
					}
				}
			}
		}
	}

	return std::make_pair(points, distances);
}