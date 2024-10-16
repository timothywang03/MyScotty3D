
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge

	FaceRef f_boundary;

	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;

	VertexRef v1 = h->next->vertex;
	VertexRef v2 = t->next->vertex;
	int v1_deg = v1->degree();
	int v2_deg = v2->degree();

	FaceRef f1 = h->face;
	FaceRef f2 = t->face;

	// create a new vertex and connect it to v1, v2
	VertexRef mid = emplace_vertex();
	EdgeRef a = emplace_edge();
	EdgeRef b = emplace_edge();
	a->sharp = e->sharp;
	b->sharp = e->sharp;
	HalfedgeRef a1 = emplace_halfedge();
	HalfedgeRef a2 = emplace_halfedge();
	HalfedgeRef b1 = emplace_halfedge();
	HalfedgeRef b2 = emplace_halfedge();
	HalfedgeRef d1, d2;

	mid->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, mid);
	interpolate_data({h, h->next}, a1);
	interpolate_data({t, t->next}, b2);
	mid->halfedge = a1;

	// create new faces
	FaceRef fa = emplace_face();
	FaceRef fb;

	// instantiate edges
	a->halfedge = a1;
	b->halfedge = b1;	

	// create new edges that will connect to mid ccw
	EdgeRef c = emplace_edge();
	EdgeRef d;
	HalfedgeRef c1 = emplace_halfedge();
	HalfedgeRef c2 = emplace_halfedge();
	c->halfedge = c1;


	// connect halfedges to new vertex
	a1->set_tnvef(a2, h->next, mid, a, fa);
	if (e->on_boundary()) {
		if (e->halfedge->face->boundary) {
			f_boundary = e->halfedge->face;
		} else {
			f_boundary = e->halfedge->twin->face;
		}
		a2->set_tnvef(a1, b2, v1, a, f2);
		b2->set_tnvef(b1, t->next, mid, b, f_boundary);
	} else {
		fb = emplace_face();
		d = emplace_edge();
		d1 = emplace_halfedge();
		d2 = emplace_halfedge();
		d1->set_tnvef(d2, t->next->next, mid, d, f2);
		d2->set_tnvef(d1, b2, t->next->next->vertex, d, fb);
		d->halfedge = d1;
		a2->set_tnvef(a1, d1, v1, a, f_boundary);
		b2->set_tnvef(b1, t->next, mid, b, fb);
	}
	b1->set_tnvef(b2, c1, v2, b, f1);
	c1->set_tnvef(c2, h->next->next, mid, c, f1);
	c2->set_tnvef(c1, a1, h->next->next->vertex, c, fa);

	v1->halfedge = a2;
	v2->halfedge = b1;

	// adjust next halfedges for the adjacents of new edges
	h->next->next = c2;
	HalfedgeRef temp = h;
	for (int i = 0; i < v1_deg - 1; i++) {
		temp = temp->next->twin;
	}
	temp->next = a2;

	temp = t;
	for (int i = 0; i < v2_deg - 1; i++) {
		temp = temp->next->twin;
	}
	temp->next = b1;

	if (!e->on_boundary()) {
		t->next->next = d2;
	}

	// iterate through to update faces
	temp = b1;
	do {
		temp->face = f1;
		temp = temp->next;
	} while (temp != b1);
	f1->halfedge = b1;

	temp = a1;
	do {
		temp->face = fa;
		temp = temp->next;
	} while (temp != a1);
	fa->halfedge = a1;

	if (!e->on_boundary()){
		temp = a2;
		do {
            temp->face = f2;
            temp = temp->next;
		} while (temp != a2);
		f2->halfedge = a2;

		temp = b2;
		do {
			temp->face = fb;
			temp = temp->next;
		} while (temp != b2);
		fb->halfedge = b2;
	} else {
		f_boundary->halfedge = a2;
	}

	// remove old halfedges
	erase_halfedge(h);
	erase_halfedge(t);
	erase_edge(e);

    return mid;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	// create a new face that is a duplicate version of the input face
	// FaceRef f_new = emplace_face();
	FaceRef f_boundary = f->halfedge->face;
	if (!f_boundary->boundary) {
		f_boundary = f->halfedge->twin->face;
	}

	HalfedgeRef h = f->halfedge;
	HalfedgeRef t = h->twin;

	std::cout << h->id << " " << t->id << std::endl;

	// create new points for the face
	std::vector<VertexRef> old_vertices;
	std::vector<VertexRef> new_vertices;
	HalfedgeRef temp = h;
	do {
		VertexRef v = emplace_vertex();
		v->position = temp->vertex->position;
		old_vertices.push_back(temp->vertex);
		new_vertices.push_back(v);
		temp = temp->next;
	} while (temp != h);

	// create new halfedges and edges for the face itself
	for (size_t i = 0; i < new_vertices.size(); i++) {
		HalfedgeRef h_i = emplace_halfedge();
		HalfedgeRef t_i = emplace_halfedge();
		EdgeRef e = emplace_edge();
		e->halfedge = h_i;
		h_i->edge = e;
		t_i->edge = e;
		h_i->twin = t_i;
		t_i->twin = h_i;
		h_i->vertex = new_vertices[i];
		t_i->vertex = new_vertices[(i + 1) % new_vertices.size()];
		new_vertices[i]->halfedge = h_i;
		h_i->face = f;
		t_i->face = f_boundary;
	}

	std::vector<HalfedgeRef> outer_edges;
	// connect the new halfedges to each other
	for (size_t i = 0; i < new_vertices.size(); i++) {
		HalfedgeRef h_i = new_vertices[i]->halfedge;
		h_i->next = new_vertices[(i + 1) % new_vertices.size()]->halfedge;
		outer_edges.push_back(h_i->twin);
	}

	f->halfedge = new_vertices[0]->halfedge;

	std::vector<HalfedgeRef> outer_joints;
	std::vector<HalfedgeRef> inner_joints;

	// connect all the vertices from the new face to the old face with halfedges
	for (size_t i = 0; i < new_vertices.size(); i++) {
		EdgeRef e = emplace_edge();
		HalfedgeRef h_i = emplace_halfedge();
		HalfedgeRef t_i = emplace_halfedge();
		e->halfedge = h_i;
		h_i->edge = e;
		t_i->edge = e;
		h_i->twin = t_i;
		t_i->twin = h_i;
		h_i->vertex = old_vertices[i];
		t_i->vertex = new_vertices[i];
		h_i->face = f_boundary;
		t_i->face = f_boundary;
		outer_joints.push_back(h_i);
		inner_joints.push_back(t_i);
	}

	// set next pointers for the new halfedges
	for (size_t i = 0; i < new_vertices.size(); i++) {
		outer_joints[i]->next = outer_edges[(i - 1) % new_vertices.size()];
		outer_edges[i]->next = inner_joints[i];
		inner_joints[i]->next = old_vertices[i]->halfedge;
		old_vertices[i]->halfedge->next = outer_joints[(i + 1) % new_vertices.size()];
	}

	// create new faces with the newly connected halfedges
	for (size_t i = 0; i < new_vertices.size(); i++) {
		FaceRef f_ = emplace_face();
		f_->halfedge = outer_joints[i];
		inner_joints[(i - 1) % new_vertices.size()]->face = f_;
		old_vertices[(i - 1) % new_vertices.size()]->halfedge->face = f_;
		outer_edges[(i - 1) % new_vertices.size()]->face = f_;
		outer_joints[i]->face = f_;
	}

	return f;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
	//A2L1: Flip Edge

	if (e->on_boundary()) return std::nullopt;

	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;

	// define vertices that are affected by the flip
	VertexRef v1 = h->next->vertex;
	VertexRef v2 = t->next->vertex;
	VertexRef v3 = h->next->next->vertex;
	VertexRef v4 = t->next->next->vertex;

	// if v1 or v2 are degree 2, return nullopt
	if (v1->degree() == 2 || v2->degree() == 2) return std::nullopt;

	// check if the flip would add duplicate halfedge
	if (h->next->next == t->next->next->twin) return std::nullopt;

	// define new faces with the new halfedges
	FaceRef f1 = h->face;
	FaceRef f2 = t->face;

	// perform the flip operation
	v1->halfedge = h->next;
	v2->halfedge = t->next;
	f1->halfedge = h;
	f2->halfedge = t;

	t->vertex = v3;
	h->vertex = v4;

	HalfedgeRef temp = h;
	for (uint32_t i = 0; i < h->face->degree() - 1; i++) {
		temp = temp->next;
	}
	temp->next = t->next;
	t->next->face = f1;
	temp = t->next->next;
	t->next->next = h;
	t->next = temp;

	temp = t;
	for (uint32_t i = 0; i < t->face->degree() - 1; i++) {
		temp = temp->next;
	}
	temp->next = h->next;
	h->next->face = f2;
	temp = h->next->next;
	h->next->next = t;
	h->next = temp;
	
	return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

	// HalfedgeRef h = v->halfedge;
	// FaceRef f_new = emplace_face();

	// HalfedgeRef temp = h->twin->next;
	// HalfedgeRef last = h;
	// std::vector<HalfedgeRef> removals;

	// do {
	// 	removals.push_back(temp);
	// 	temp->next->next = last->next;
	// 	erase_face(temp->face);
	// 	last = temp;
	// } while (temp != h->twin->next);

	// temp->next->next = last->next;

	// for (size_t i = 0; i < removals.size(); i++) {
	// 	erase_halfedge(removals[i]->twin);
	// 	erase_halfedge(removals[i]);
	// 	erase_edge(removals[i]->edge);
	// }

	// temp = h->next;
	// do {
	// 	temp->face = f_new;
	// 	temp = temp->next;
	// } while (temp != h->next);

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

	// if any of the vertices neighboring any of the collapsed edge vertices have degree 2, then return nullopt

	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;

	VertexRef v1 = h->next->vertex;
	VertexRef v2 = t->next->vertex;
	int v1_degree = v1->degree();
	int v2_degree = v2->degree();

	// creates a new vertex in between v1 and v2
	VertexRef x = emplace_vertex();
	x->position = (v1->position + v2->position) / 2.0f;
	x->halfedge = h->next;
	interpolate_data({v1, v2}, x);

	// change halfedge references
	h->next->vertex = x;
	t->next->vertex = x;	

	bool hflag = false, tflag = false;

	HalfedgeRef temp = h->next->next->twin->next;
	HalfedgeRef redundant = h->next->next;
	if (h->next->next->next == h) {
		h->next->twin->vertex->halfedge = h->next->twin;
		erase_edge(redundant->edge);
		erase_halfedge(redundant->twin);
		erase_halfedge(redundant);
		h->next->next = temp;

		temp = t;
		for (uint32_t i = 0; i < h->next->vertex->degree() - 2; i++) {
			temp = temp->next->twin;
		}
		temp->next = h->next;
		erase_face(h->face);
		erase_face(h->next->next->face);
		
		// delete two faces and iterate through to make new one
		FaceRef joined_face = emplace_face();
		temp = h->next;
		do {
			temp->face = joined_face;
			temp = temp->next;
		} while (temp != h->next);
		joined_face->halfedge = h->next;
		hflag = true;
	}

	temp = t->next->next->twin->next;
	redundant = h->next->next;
	if (t->next->next->next == t) {
		t->next->twin->vertex->halfedge = t->next->twin;
		erase_edge(redundant->edge);
		erase_halfedge(redundant->twin);
		erase_halfedge(redundant);
		t->next->next = temp;

		temp = h;
		for (uint32_t i = 0; i < t->next->vertex->degree() - 2; i++) {
			temp = temp->next->twin;
		}
		temp->next = t->next;
		erase_face(t->face);
		erase_face(t->next->next->face);

		// delete two faces and iterate through to make new one
		FaceRef joined_face = emplace_face();
		temp = t->next;
		do {
			temp->face = joined_face;
			temp = temp->next;
		} while (temp != t->next);
		joined_face->halfedge = t->next;
		tflag = true;
	}

	if (!tflag) {
		temp = h;
		for (int i = 0; i < v1_degree - 1; i++) {
			temp = temp->next->twin;
			temp->twin->vertex = x;
		}
		temp->next = t->next;
	}

	if (!hflag) {
		temp = t;
		for (int i = 0; i < v2_degree - 1; i++) {
			temp = temp->next->twin;
			temp->twin->vertex = x;
		}
		temp->next = h->next;
	}

	h->face->halfedge = h->next;
	t->face->halfedge = t->next;

	erase_vertex(v1);
	erase_vertex(v2);
	erase_halfedge(h);
	erase_halfedge(t);
	erase_edge(e);
	
    return x;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move

	HalfedgeRef h = face->halfedge;
	Vec3 center = face->center();
	HalfedgeRef temp = h;
	do {
		temp->vertex->position = temp->vertex->position * (1.0f - shrink) + center * shrink;
		temp->vertex->position = temp->vertex->position + move;
		temp = temp->next;
	} while (temp != h);
}
