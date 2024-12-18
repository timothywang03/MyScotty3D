#include "test.h"
#include "geometry/halfedge.h"
#include <iostream>

static void expect_split(Halfedge_Mesh &mesh, Halfedge_Mesh::EdgeRef edge, Halfedge_Mesh const &after) {
	if (auto ret = mesh.split_edge(edge)) {
		if (auto msg = mesh.validate()) {
			throw Test::error("Invalid mesh: " + msg.value().second);
		}
		// check mesh shape:
		if (auto difference = Test::differs(mesh, after, Test::CheckAllBits)) {
			throw Test::error("Resulting mesh did not match expected: " + *difference);
		}
	} else {
		throw Test::error("split_edge rejected operation!");
	}
}

/*
BASIC CASE:

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Split Edge on Edge: 1-4

After mesh:
0--1\
|\ | \
| \2--3
|  | /
4--5/
*/
Test test_a2_l2_split_edge_basic_simple("a2.l2.split_edge.basic.simple", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                         Vec3(1.25f, 0.0f, 0.0f),  Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 4, 5, 2}, 
		{0, 2, 1}, 
		{1, 2, 3}, 
		{2, 5, 3}
	});

	expect_split(mesh, edge, after);
});

/*
EDGE CASE: 

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Split Edge on Edge: 0-1

After mesh:
0--1--2\
|  /  | \
| /   |  3
|/    | /
4-----5/
*/
Test test_a2_l2_split_edge_edge_boundary("a2.l2.split_edge.edge.boundary", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f),  Vec3(0.05f, 1.05f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            						Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), 							Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 4, 1}, 
		{1, 4, 5, 2}, 
		{2, 5, 3}
	});

	expect_split(mesh, edge, after);
});

/*
Edge CASE

Initial mesh:
0--1--2
|  |  |
|  |  |
|  |  |
|  |  |
|  |  |
3--4--5

Split Edge on Edge: 1-4

After mesh:
0--1--2
|\ |  |
| \|  |
|  6  |
|  |\ |
|  | \|
3--4--5
*/
Test test_a2_l2_split_edge_edge_square_diagonal(
    "a2.l2.split_edge.square.diagonal", []() {
      Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
          {Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 0.0f, 0.0f),
           Vec3(4.0f, 0.0f, 0.0f), Vec3(0.0f, 4.0f, 0.0f),
           Vec3(2.0f, 4.0f, 0.0f), Vec3(4.0f, 4.0f, 0.0f)},
          {{0, 3, 4, 1}, {1, 4, 5, 2}});
      Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

      Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces(
          {Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 0.0f, 0.0f),
           Vec3(4.0f, 0.0f, 0.0f), Vec3(0.0f, 4.0f, 0.0f),
           Vec3(2.0f, 4.0f, 0.0f), Vec3(4.0f, 4.0f, 0.0f),
           Vec3(2.0f, 2.0f, 0.0f)},
          {{1, 0, 6}, {1, 6, 5, 2}, {4, 5, 6}, {0, 3, 4, 6}});

      expect_split(mesh, edge, after);
    });

// triangle case
/*
Edge CASE

Initial mesh:
1
|\
| \
|  \
0---2

Split Edge on Edge: 1-2

After mesh:
1
|\
| \
|  3
| / \
|/   \
0-----2
*/

Test test_a2_l2_split_edge_edge_triangle("a2.l2.split_edge.triangle", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3 (0.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f)}, {{1, 0, 2}});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3 (0.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(0.5f, 0.5f, 0.0f)}, 
		{{1, 0, 3}, {3, 0, 2}});

	expect_split(mesh, edge, mesh);
});

/*
Edge CASE

Initial mesh:
1
|\
| \
|  \
0---2

Split Edge on Edge: 1-2

After mesh:
1
|\
| \
|\ \
| \ \
|  \ \
0--3--2
*/

Test test_a2_l2_split_edge_edge_triangle2("a2.l2.split_edge.triangle2", []() {
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f)},
      {{1, 0, 2}});
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

  Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces(
      {Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
       Vec3(0.5f, 0.5f, 0.0f)},
      {{1, 0, 3}, {1, 3, 2}});

  expect_split(mesh, edge, mesh);
});

/*
Edge CASE

Initial mesh:
0--1
|  |
|  |
|  |
|  |
|  |
2--3

Split Edge on Edge: 1-3

After mesh:
0--1
|\ |
| \|
|  4
|  |
|  |
2--3
*/
Test test_a2_l2_split_edge_edge_square("a2.l2.split_edge.square", []() {
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 0.0f, 0.0f), Vec3(0.0f, 4.0f, 0.0f),
       Vec3(2.0f, 4.0f, 0.0f)},
      {{0, 2, 3, 1}});
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

  Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces(
      {Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 0.0f, 0.0f), Vec3(0.0f, 4.0f, 0.0f),
       Vec3(2.0f, 4.0f, 0.0f), Vec3(2.0f, 2.0f, 0.0f)},
      {{0, 4, 1}, {0, 2, 3, 4}});

  expect_split(mesh, edge, after);
});

Test test1("a2.l2.split_edge.isosceles_triangle", []() {
  // Create an isosceles triangle mesh
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {
          Vec3(0.0f, 0.0f, 0.0f),  // Vertex 0: left base corner
          Vec3(1.0f, 0.0f, 0.0f),  // Vertex 1: right base corner
          Vec3(0.5f, std::sqrt(3.0f) / 2.0f, 0.0f)  // Vertex 2: apex
      },
      {{0, 1, 2}});

  // Get the hypotenuse edge (the edge between vertices 1 and 2)
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

  // Create the expected mesh after splitting
  Halfedge_Mesh expected = Halfedge_Mesh::from_indexed_faces(
      {
          Vec3(0.0f, 0.0f, 0.0f),  // Vertex 0: left base corner
          Vec3(1.0f, 0.0f, 0.0f),  // Vertex 1: right base corner
          Vec3(0.5f, std::sqrt(3.0f) / 2.0f, 0.0f),  // Vertex 2: apex
          Vec3(0.75f, std::sqrt(3.0f) / 4.0f,
               0.0f)  // Vertex 3: new vertex on hypotenuse
      },
      {{0, 1, 3}, {0, 3, 2}});
  expect_split(mesh, edge, expected);
});

Test test2("a2.l2.split_edge.right_triangle", []() {
  // Create a right triangle mesh (3-4-5 triangle)
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {
          Vec3(0.0f, 0.0f, 0.0f),  // Vertex 0: bottom-left corner
          Vec3(4.0f, 0.0f, 0.0f),  // Vertex 1: bottom-right corner
          Vec3(0.0f, 3.0f, 0.0f)   // Vertex 2: top-left corner
      },
      {{0, 1, 2}});

  // Get the hypotenuse edge (the edge between vertices 1 and 2)
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

  // Create the expected mesh after splitting
  Halfedge_Mesh expected = Halfedge_Mesh::from_indexed_faces(
      {
          Vec3(0.0f, 0.0f, 0.0f),  // Vertex 0: bottom-left corner
          Vec3(4.0f, 0.0f, 0.0f),  // Vertex 1: bottom-right corner
          Vec3(0.0f, 3.0f, 0.0f),  // Vertex 2: top-left corner
          Vec3(2.0f, 1.5f, 0.0f)   // Vertex 3: new vertex on hypotenuse
      },
      {{0, 1, 3}, {0, 3, 2}});

  expect_split(mesh, edge, expected);

});