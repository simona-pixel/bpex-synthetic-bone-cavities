# Digital hollow shaft component

The design task is a hollow shaft component for a modular training phantom. The source is the outer envelope of public Figshare tibia 102480, Keast, Bonacci and Fox, Statistical Shape Model of the Tibia version 4, doi:10.6084/m9.figshare.20454462.v4, CC BY 4.0. See ATTRIBUTION.md for source and derivative attribution.

The preserved v3 profile constructs and checks a synthetic cavity. A central 40-mm box clips the outer geometry. Subtraction of the longer cavity opens a through lumen with annular end faces. The supplied figshare_102480_hollow_40mm.stl is the resulting material boundary, in millimetres. This derivative supplies a digital design example, not an imaged endosteal surface.

Observed digital checks: one connected component, zero boundary or non-manifold edges, consistent winding, no detected non-adjacent self-intersections, and Euler characteristic zero. Material and removed volumes were 14.259 and 3.979 cubic centimetres. Full numerical values are in ../results/qualification/digital_coupon_report.json. Parent-surface clearance is not a wall-thickness guarantee at the clipped rims.

The component can be inspected in a mesh viewer without a printer. Insert dimensions, assembly, material, printing, instrumentation, mechanical and radiographic performance were not evaluated.
