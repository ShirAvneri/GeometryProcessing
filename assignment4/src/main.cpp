#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>

#include <igl/cotmatrix.h>
#include <igl/slice.h>

#include "Lasso.h"
#include "Colors.h"

//activate this for alternate UI (easier to debug)
#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

//vertex array, #V x3
Eigen::MatrixXd V(0, 3), V_cp(0, 3);
//face array, #F x3
Eigen::MatrixXi F(0, 3);

//mouse interaction
enum MouseMode { SELECT, TRANSLATE, ROTATE, NONE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0, 1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0, 1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0, 1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0, 3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0, 3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0, 0, 0);
Eigen::Vector4f rotation(0, 0, 0, 1.);
typedef Eigen::Triplet<double> T;
//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

//function declarations (see below for implementation)
bool solve(Viewer& viewer);
void get_new_handle_locations();
Eigen::Vector3f computeTranslation(Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
Eigen::Vector4f computeRotation(Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
void compute_handle_centroids();
Eigen::MatrixXd readMatrix(const char* filename);

bool callback_mouse_down(Viewer& viewer, int button, int modifier);
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(Viewer& viewer, int button, int modifier);
bool callback_pre_draw(Viewer& viewer);
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void onNewHandleID();
void applySelection();

typedef SparseMatrix<double> SparseXd;
enum BREAKPOINT { S = 0, B = 1, B_TAG = 2, S_TAG = 3 };
int breakpoint = BREAKPOINT::S;
string mesh_path = "";

bool solve(Viewer& viewer)
{
	/**** Add your code for computing the deformation from handle_vertex_positions and handle_vertices here (replace following line) ****/
	static int prev_handle = -1;
	static SimplicialCholesky<SparseXd, RowMajor> solver;
	static SparseXd A, A_ff, A_fc;

	size_t vertices_count = V.rows(), handles_count = handle_vertices.rows();

	VectorXi free_vertices(vertices_count - handles_count);
	size_t i = 0;
	for (size_t j = 0; j < vertices_count; j++)
	{
		if (handle_id(j) == -1)
		{
			free_vertices(i) = j;
			i++;
		}
	}

	/*** STEP 2 ***/
	if (moving_handle != prev_handle)
	{ // Need to update A, A_ff, A_fc
		SparseXd L, M;
		igl::cotmatrix(V, F, L);
		igl::massmatrix(V, F, igl::MassMatrixType(), M);

		A = L * M.cwiseInverse() * L;

		igl::slice(A, free_vertices, free_vertices, A_ff);
		igl::slice(A, free_vertices, handle_vertices, A_fc);

		solver.compute(A_ff);

		prev_handle = moving_handle;
	}

	if (breakpoint == BREAKPOINT::S)
		return true;

	MatrixXd original_positions(handles_count, 1);
	igl::slice(V, handle_vertices, 1, original_positions);
	igl::slice_into(solver.solve(-A_fc * original_positions), free_vertices, 1, V);

	// V_S --> V_cp
	// V_B --> V



	/*** STEP 3 ***/
	if (breakpoint == BREAKPOINT::B)
		return true;

	static MatrixXd d(V_cp - V);

	MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	vector<vector<int>> VV;
	igl::adjacency_list(F, VV);

	MatrixXd x(vertices_count, 3);
	VectorXi longest_edges(vertices_count);
	for (size_t i = 0; i < vertices_count; i++)
	{
		RowVector3d longest_proj_e(0, 0, 0);
		size_t neighbors_count = VV[i].size();
		for (size_t j = 0; j < neighbors_count; j++)
		{
			RowVector3d e = V.row(VV[i][j]) - V.row(i);
			RowVector3d proj_e = e - e.dot(N.row(i)) * N.row(i);
			if (longest_proj_e.norm() < proj_e.norm())
			{
				longest_proj_e = proj_e;
				longest_edges(i) = VV[i][j];
			}
		}

		x.row(i) = longest_proj_e.normalized();
	}

	MatrixXd y(vertices_count, 3);
	for (size_t i = 0; i < vertices_count; i++)
	{
		y.row(i) = RowVector3d(N.row(i)).cross(RowVector3d(x.row(i)));
	}

	for (size_t i = 0; i < vertices_count; i++) {
		RowVector3d d_i(d.row(i)), x_i(x.row(i)), y_i(y.row(i)), N_i(N.row(i));
		d.row(i) = d_i.dot(x_i) * x_i + d_i.dot(y_i) * y_i + d_i.dot(N_i) * N_i;
	}

	/*** STEP 4 ***/
	igl::slice_into(solver.solve(-A_fc * handle_vertex_positions), free_vertices, 1, V);
	igl::slice_into(handle_vertex_positions, handle_vertices, 1, V);

	if (breakpoint == BREAKPOINT::B_TAG)
		return true;

	/*** STEP 5 ***/
	igl::per_vertex_normals(V, F, N);

	for (size_t i = 0; i < vertices_count; i++)
	{
		RowVector3d e = V.row(longest_edges(i)) - V.row(i);
		x.row(i) = (e - e.dot(N.row(i)) * N.row(i)).normalized();
	}

	for (size_t i = 0; i < vertices_count; i++)
	{
		y.row(i) = RowVector3d(N.row(i)).cross(RowVector3d(x.row(i)));
	}

	static MatrixXd transformed_d(vertices_count, 3);
	for (size_t i = 0; i < vertices_count; i++) {
		RowVector3d d_i(d.row(i)), x_i(x.row(i)), y_i(y.row(i)), N_i(N.row(i));
		transformed_d.row(i) = d_i.dot(x_i) * x_i + d_i.dot(y_i) * y_i + d_i.dot(N_i) * N_i;
	}

	V += transformed_d;
	V_cp = V;

	/*** DONE! ***/
	return true;
};

void get_new_handle_locations()
{
	int count = 0;
	for (long vi = 0; vi < V.rows(); ++vi)
		if (handle_id[vi] >= 0)
		{
			Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();

			if (handle_id[vi] == moving_handle) {
				if (mouse_mode == TRANSLATE)
					goalPosition += translation;
				else if (mouse_mode == ROTATE) {
					Eigen::RowVector3f  goalPositionCopy = goalPosition;
					goalPosition -= handle_centroids.row(moving_handle).cast<float>();
					igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPositionCopy.data());
					goalPosition = goalPositionCopy;
					goalPosition += handle_centroids.row(moving_handle).cast<float>();
				}
			}
			handle_vertex_positions.row(count++) = goalPosition.cast<double>();
		}
}

bool load_mesh(string filename)
{
	mesh_path = filename;

	igl::read_triangle_mesh(filename, V, F);
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	viewer.core.align_camera_center(V);
	V_cp = V;
	handle_id.setConstant(V.rows(), 1, -1);
	// Initialize selector
	lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

	selected_v.resize(0, 1);

	return true;
}

int main(int argc, char* argv[])
{
	if (argc != 2) 
	{
		cout << "Usage assignment4 mesh.off>" << endl;
		load_mesh("../data/camel_head.off");
	}
	else
	{
		load_mesh(argv[1]);
	}

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Deformation Controls", ImGuiTreeNodeFlags_DefaultOpen))
		{
			int mouse_mode_type = static_cast<int>(mouse_mode);

			if (ImGui::Combo("Mouse Mode", &mouse_mode_type, "SELECT\0TRANSLATE\0ROTATE\0NONE\0"))
			{
				mouse_mode = static_cast<MouseMode>(mouse_mode_type);
			}

			if (ImGui::Button("Clear Selection", ImVec2(-1, 0)))
			{
				selected_v.resize(0, 1);
			}

			if (ImGui::Button("Apply Selection", ImVec2(-1, 0)))
			{
				applySelection();
			}

			if (ImGui::Button("Clear Constraints", ImVec2(-1, 0)))
			{
				handle_id.setConstant(V.rows(), 1, -1);
			}

			if (ImGui::Button("Reset", ImVec2(-1, 0)))
			{
				load_mesh(mesh_path);
			}

			if (ImGui::TreeNode("Breakpoints"))
			{
				if (ImGui::RadioButton("S", &breakpoint, (int)BREAKPOINT::S) ||
					ImGui::RadioButton("B", &breakpoint, (int)BREAKPOINT::B) ||
					ImGui::RadioButton("B'", &breakpoint, (int)BREAKPOINT::B_TAG) ||
					ImGui::RadioButton("S'", &breakpoint, (int)BREAKPOINT::S_TAG))
				{
					solve(viewer);
				}
				ImGui::TreePop();
			}
		}
	};

	viewer.callback_key_down = callback_key_down;
	viewer.callback_mouse_down = callback_mouse_down;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_mouse_up = callback_mouse_up;
	viewer.callback_pre_draw = callback_pre_draw;

	viewer.data().point_size = 10;
	viewer.core.set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
	viewer.launch();
}

bool callback_mouse_down(Viewer& viewer, int button, int modifier)
{
	if (button == (int)Viewer::MouseButton::Right)
		return false;

	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;

	if (mouse_mode == SELECT)
	{
		if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >= 0)
			doit = true;
		else
			lasso->strokeReset();
	}
	else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
	{
		int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
		if (vi >= 0 && handle_id[vi] >= 0)  //if a region was found, mark it for translation/rotation
		{
			moving_handle = handle_id[vi];
			get_new_handle_locations();
			doit = true;
		}
	}
	return doit;
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y)
{
	if (!doit)
		return false;
	if (mouse_mode == SELECT)
	{
		lasso->strokeAdd(mouse_x, mouse_y);
		return true;
	}
	if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
	{
		if (mouse_mode == TRANSLATE) {
			translation = computeTranslation(viewer,
				mouse_x,
				down_mouse_x,
				mouse_y,
				down_mouse_y,
				handle_centroids.row(moving_handle));
		}
		else {
			rotation = computeRotation(viewer,
				mouse_x,
				down_mouse_x,
				mouse_y,
				down_mouse_y,
				handle_centroids.row(moving_handle));
		}
		get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
		solve(viewer);
		down_mouse_x = mouse_x;
		down_mouse_y = mouse_y;
#endif
		return true;

		}
	return false;
	}

bool callback_mouse_up(Viewer& viewer, int button, int modifier)
{
	if (!doit)
		return false;
	doit = false;
	if (mouse_mode == SELECT)
	{
		selected_v.resize(0, 1);
		lasso->strokeFinish(selected_v);
		return true;
	}

	if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
	{
#ifdef UPDATE_ONLY_ON_UP
		if (moving_handle >= 0)
			solve(viewer);
#endif
		translation.setZero();
		rotation.setZero(); rotation[3] = 1.;
		moving_handle = -1;

		compute_handle_centroids();

		return true;
	}

	return false;
};


bool callback_pre_draw(Viewer& viewer)
{
	// initialize vertex colors
	vertex_colors = Eigen::MatrixXd::Constant(V.rows(), 3, .9);

	// first, color constraints
	int num = handle_id.maxCoeff();
	if (num == 0)
		num = 1;
	for (int i = 0; i < V.rows(); ++i)
		if (handle_id[i] != -1)
		{
			int r = handle_id[i] % MAXNUMREGIONS;
			vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
		}
	// then, color selection
	for (int i = 0; i < selected_v.size(); ++i)
		vertex_colors.row(selected_v[i]) << 131. / 255, 131. / 255, 131. / 255.;

	viewer.data().set_colors(vertex_colors);
	viewer.data().V_material_specular.fill(0);
	viewer.data().V_material_specular.col(3).fill(1);
	viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;


	//clear points and lines
	viewer.data().set_points(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));
	viewer.data().set_edges(Eigen::MatrixXd::Zero(0, 3), Eigen::MatrixXi::Zero(0, 3), Eigen::MatrixXd::Zero(0, 3));

	//draw the stroke of the selection
	for (unsigned int i = 0; i < lasso->strokePoints.size(); ++i)
	{
		viewer.data().add_points(lasso->strokePoints[i], Eigen::RowVector3d(0.4, 0.4, 0.4));
		if (i > 1)
			viewer.data().add_edges(lasso->strokePoints[i - 1], lasso->strokePoints[i], Eigen::RowVector3d(0.7, 0.7, 0.7));
	}

	// update the vertex position all the time
	viewer.data().V.resize(V.rows(), 3);
	viewer.data().V << V;

	viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_POSITION;

#ifdef UPDATE_ONLY_ON_UP
	//draw only the moving parts with a white line
	if (moving_handle >= 0)
	{
		Eigen::MatrixXd edges(3 * F.rows(), 6);
		int num_edges = 0;
		for (int fi = 0; fi < F.rows(); ++fi)
		{
			int firstPickedVertex = -1;
			for (int vi = 0; vi < 3; ++vi)
				if (handle_id[F(fi, vi)] == moving_handle)
				{
					firstPickedVertex = vi;
					break;
				}
			if (firstPickedVertex == -1)
				continue;


			Eigen::Matrix3d points;
			for (int vi = 0; vi < 3; ++vi)
			{
				int vertex_id = F(fi, vi);
				if (handle_id[vertex_id] == moving_handle)
				{
					int index = -1;
					// if face is already constrained, find index in the constraints
					(handle_vertices.array() - vertex_id).cwiseAbs().minCoeff(&index);
					points.row(vi) = handle_vertex_positions.row(index);
				}
				else
					points.row(vi) = V.row(vertex_id);

			}
			edges.row(num_edges++) << points.row(0), points.row(1);
			edges.row(num_edges++) << points.row(1), points.row(2);
			edges.row(num_edges++) << points.row(2), points.row(0);
		}
		edges.conservativeResize(num_edges, Eigen::NoChange);
		viewer.data().add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9, 0.9, 0.9));

	}
#endif
	return false;

}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers)
{
	bool handled = false;
	if (key == 'S')
	{
		mouse_mode = SELECT;
		handled = true;
	}

	if ((key == 'T') && (modifiers == IGL_MOD_ALT))
	{
		mouse_mode = TRANSLATE;
		handled = true;
	}

	if ((key == 'R') && (modifiers == IGL_MOD_ALT))
	{
		mouse_mode = ROTATE;
		handled = true;
	}
	if (key == 'A')
	{
		applySelection();
		callback_key_down(viewer, '1', 0);
		handled = true;
	}
	// if (key == '2')
	// viewer.core.align_camera_position(B);
	//  }

	//viewer.ngui->refresh();
	return handled;
}

void onNewHandleID()
{
	//store handle vertices too
	int numFree = (handle_id.array() == -1).cast<int>().sum();
	int num_handle_vertices = V.rows() - numFree;
	handle_vertices.setZero(num_handle_vertices);
	handle_vertex_positions.setZero(num_handle_vertices, 3);

	int count = 0;
	for (long vi = 0; vi < V.rows(); ++vi)
		if (handle_id[vi] >= 0)
			handle_vertices[count++] = vi;

	compute_handle_centroids();
}

void applySelection()
{
	int index = handle_id.maxCoeff() + 1;
	for (int i = 0; i < selected_v.rows(); ++i)
	{
		const int selected_vertex = selected_v[i];
		if (handle_id[selected_vertex] == -1)
			handle_id[selected_vertex] = index;
	}
	selected_v.resize(0, 1);

	onNewHandleID();
}

void compute_handle_centroids()
{
	//compute centroids of handles
	int num_handles = handle_id.maxCoeff() + 1;
	handle_centroids.setZero(num_handles, 3);

	Eigen::VectorXi num; num.setZero(num_handles, 1);
	for (long vi = 0; vi < V.rows(); ++vi)
	{
		int r = handle_id[vi];
		if (r != -1)
		{
			handle_centroids.row(r) += V.row(vi);
			num[r]++;
		}
	}

	for (long i = 0; i < num_handles; ++i)
		handle_centroids.row(i) = handle_centroids.row(i).array() / num[i];

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f computeTranslation(Viewer& viewer,
	int mouse_x,
	int from_x,
	int mouse_y,
	int from_y,
	Eigen::RowVector3d pt3D)
{
	Eigen::Matrix4f modelview = viewer.core.view;// * viewer.data().model;
	//project the given point (typically the handle centroid) to get a screen space depth
	Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
		modelview,
		viewer.core.proj,
		viewer.core.viewport);
	float depth = proj[2];

	double x, y;
	Eigen::Vector3f pos1, pos0;

	//unproject from- and to- points
	x = mouse_x;
	y = viewer.core.viewport(3) - mouse_y;
	pos1 = igl::unproject(Eigen::Vector3f(x, y, depth),
		modelview,
		viewer.core.proj,
		viewer.core.viewport);


	x = from_x;
	y = viewer.core.viewport(3) - from_y;
	pos0 = igl::unproject(Eigen::Vector3f(x, y, depth),
		modelview,
		viewer.core.proj,
		viewer.core.viewport);

	//translation is the vector connecting the two
	Eigen::Vector3f translation = pos1 - pos0;
	return translation;

}


//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f computeRotation(Viewer& viewer,
	int mouse_x,
	int from_x,
	int mouse_y,
	int from_y,
	Eigen::RowVector3d pt3D)
{

	Eigen::Vector4f rotation;
	rotation.setZero();
	rotation[3] = 1.;

	Eigen::Matrix4f modelview = viewer.core.view;// * viewer.data().model;

	//initialize a trackball around the handle that is being rotated
	//the trackball has (approximately) width w and height h
	double w = viewer.core.viewport[2] / 8;
	double h = viewer.core.viewport[3] / 8;

	//the mouse motion has to be expressed with respect to its center of mass
	//(i.e. it should approximately fall inside the region of the trackball)

	//project the given point on the handle(centroid)
	Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
		modelview,
		viewer.core.proj,
		viewer.core.viewport);
	proj[1] = viewer.core.viewport[3] - proj[1];

	//express the mouse points w.r.t the centroid
	from_x -= proj[0]; mouse_x -= proj[0];
	from_y -= proj[1]; mouse_y -= proj[1];

	//shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
	from_x += w / 2; mouse_x += w / 2;
	from_y += h / 2; mouse_y += h / 2;

	//get rotation from trackball
	Eigen::Vector4f drot = viewer.core.trackball_angle.coeffs();
	Eigen::Vector4f drot_conj;
	igl::quat_conjugate(drot.data(), drot_conj.data());
	igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());

	//account for the modelview rotation: prerotate by modelview (place model back to the original
	//unrotated frame), postrotate by inverse modelview
	Eigen::Vector4f out;
	igl::quat_mult(rotation.data(), drot.data(), out.data());
	igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
	return rotation;
}