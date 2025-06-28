
#!/usr/bin/env python3
"""
Interface Resultant Analyzer for Masonry Structures
Computes resultant forces and eccentricities across block interfaces
to construct lines of thrust
"""

import re
import math
import csv
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

@dataclass
class FaceForce:
    """Force data at a face vertex"""
    face_id: int
    point: Tuple[float, float, float]  # (x, y, z)
    fn: float  # Normal force
    ft: float  # Tangential force

@dataclass
class FaceGeometry:
    """Geometric data for a face"""
    face_id: int
    block_j: int
    block_k: int
    pt1: Tuple[float, float, float]
    pt2: Tuple[float, float, float]
    normal: Tuple[float, float]    # (nx, ny)
    tangent: Tuple[float, float]   # (tx, ty)
    
    @property
    def length(self) -> float:
        """Calculate face length"""
        dx = self.pt2[0] - self.pt1[0]
        dy = self.pt2[1] - self.pt1[1]
        return math.sqrt(dx*dx + dy*dy)
    
    @property
    def edge_direction(self) -> Tuple[float, float]:
        """Unit vector from pt1 to pt2"""
        dx = self.pt2[0] - self.pt1[0]
        dy = self.pt2[1] - self.pt1[1]
        length = self.length
        if length < 1e-12:
            return (0, 0)
        return (dx/length, dy/length)

@dataclass
class InterfaceResultant:
    """Resultant force and properties for an interface"""
    block_a: int
    block_b: int
    centroid: Tuple[float, float]
    virtual_face_length: float
    tangent: Tuple[float, float]      # Virtual face tangent direction
    normal: Tuple[float, float]       # Virtual face normal direction
    resultant_fx: float               # Global X force
    resultant_fy: float               # Global Y force
    resultant_magnitude: float
    moment_z: float                   # Moment about centroid
    eccentricity: float               # Distance from centroid
    application_point: Tuple[float, float]  # Where resultant acts

class InterfaceAnalyzer:
    """Main analyzer class"""
    
    def __init__(self):
        self.face_forces: List[FaceForce] = []
        self.face_geometries: Dict[int, FaceGeometry] = {}
    
    def parse_results_file(self, filepath: str) -> None:
        """Parse the results file to extract face forces and geometry"""
        print(f"Parsing results file: {filepath}")
        
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Parse face forces section
        self._parse_face_forces(content)
        
        # Parse face geometry section  
        self._parse_face_geometry(content)
        
        print(f"Loaded {len(self.face_forces)} force entries")
        print(f"Loaded {len(self.face_geometries)} face geometries")
    
    def _parse_face_forces(self, content: str) -> None:
        """Parse the === FaceID; Pt; fn; ft === section"""
        # Find the forces section
        pattern = r'=== FaceID; Pt; fn; ft ===\s*\n(.*?)\n(?=\S|\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            raise ValueError("Could not find face forces section")
        
        forces_text = match.group(1).strip()
        
        for line in forces_text.split('\n'):
            line = line.strip()
            if not line:
                continue
                
            # Parse: faceId; ptX,ptY,ptZ; fn; ft
            parts = line.split(';')
            if len(parts) != 4:
                continue
                
            face_id = int(parts[0].strip())
            
            # Parse coordinates
            coords_str = parts[1].strip()
            coords = [float(x.strip()) for x in coords_str.split(',')]
            
            fn = float(parts[2].strip())
            ft = float(parts[3].strip())
            
            self.face_forces.append(FaceForce(
                face_id=face_id,
                point=(coords[0], coords[1], coords[2]),
                fn=fn,
                ft=ft
            ))
    
    def _parse_face_geometry(self, content: str) -> None:
        """Parse the --- Face ID, BlockJ, BlockK, Pt1, Pt2, Normal, Tangent --- section"""
        # Find the geometry section
        pattern = r'--- Face ID, BlockJ, BlockK, Pt1, Pt2, Normal, Tangent ---\s*\n(.*?)\n(?=\S|\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            raise ValueError("Could not find face geometry section")
        
        geometry_text = match.group(1).strip()
        
        for line in geometry_text.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            # Parse: faceId;blockJ;blockK;x1,y1,z1; x2,y2,z2; nx,ny; tx,ty
            parts = line.split(';')
            if len(parts) != 7:
                continue
            
            face_id = int(parts[0].strip())
            block_j = int(parts[1].strip())
            block_k = int(parts[2].strip())
            
            # Parse point 1
            pt1_coords = [float(x.strip()) for x in parts[3].strip().split(',')]
            
            # Parse point 2
            pt2_coords = [float(x.strip()) for x in parts[4].strip().split(',')]
            
            # Parse normal
            normal_coords = [float(x.strip()) for x in parts[5].strip().split(',')]
            
            # Parse tangent
            tangent_coords = [float(x.strip()) for x in parts[6].strip().split(',')]
            
            self.face_geometries[face_id] = FaceGeometry(
                face_id=face_id,
                block_j=block_j,
                block_k=block_k,
                pt1=(pt1_coords[0], pt1_coords[1], pt1_coords[2]),
                pt2=(pt2_coords[0], pt2_coords[1], pt2_coords[2]),
                normal=(normal_coords[0], normal_coords[1]),
                tangent=(tangent_coords[0], tangent_coords[1])
            )
    
    def find_interfaces(self) -> List[Tuple[int, int]]:
        """Find all unique block-to-block interfaces"""
        interfaces = set()
        
        for geom in self.face_geometries.values():
            # Only consider interfaces between structural blocks (positive IDs)
            if geom.block_j > 0 and geom.block_k > 0:
                # Always order as (smaller, larger) for consistency
                block_a = min(geom.block_j, geom.block_k)
                block_b = max(geom.block_j, geom.block_k)
                interfaces.add((block_a, block_b))
        
        return sorted(list(interfaces))
    
    def analyze_interface(self, block_a: int, block_b: int) -> InterfaceResultant:
        """Analyze a specific interface between two blocks"""
        # Find all faces that connect these blocks
        interface_faces = []
        for geom in self.face_geometries.values():
            if ((geom.block_j == block_a and geom.block_k == block_b) or
                (geom.block_j == block_b and geom.block_k == block_a)):
                interface_faces.append(geom)
        
        if not interface_faces:
            raise ValueError(f"No interface found between blocks {block_a} and {block_b}")
        
        # Collect all vertices from interface faces
        all_vertices = []
        for face in interface_faces:
            all_vertices.extend([face.pt1[:2], face.pt2[:2]])  # Only X,Y coordinates
        
        # Remove duplicates and calculate centroid
        unique_vertices = list(set(all_vertices))
        centroid_x = sum(v[0] for v in unique_vertices) / len(unique_vertices)
        centroid_y = sum(v[1] for v in unique_vertices) / len(unique_vertices)
        centroid = (centroid_x, centroid_y)
        
        # Find virtual face orientation (longest edge direction)
        longest_length = 0
        virtual_tangent = (1, 0)  # Default
        
        for face in interface_faces:
            if face.length > longest_length:
                longest_length = face.length
                virtual_tangent = face.edge_direction
        
        # Virtual face normal (perpendicular to tangent)
        virtual_normal = (-virtual_tangent[1], virtual_tangent[0])
        
        # Calculate virtual face length (span of all vertices)
        min_x = min(v[0] for v in unique_vertices)
        max_x = max(v[0] for v in unique_vertices)
        min_y = min(v[1] for v in unique_vertices)
        max_y = max(v[1] for v in unique_vertices)
        virtual_length = max(max_x - min_x, max_y - min_y)
        
        # Sum forces acting on block_b
        total_fx = 0.0
        total_fy = 0.0
        total_moment = 0.0
        
        for face in interface_faces:
            # Get all forces for this face
            face_forces = [f for f in self.face_forces if f.face_id == face.face_id]
            
            for force in face_forces:
                # Convert local forces (fn, ft) to global (fx, fy)
                # Force acts on block_b, so determine sign
                if face.block_k == block_b:
                    # Normal points from J to K, so from block_a to block_b
                    sign = 1.0
                else:
                    # Normal points from block_b to block_a
                    sign = -1.0
                
                # Global force components
                fx = sign * (force.fn * face.normal[0] + force.ft * face.tangent[0])
                fy = sign * (force.fn * face.normal[1] + force.ft * face.tangent[1])
                
                total_fx += fx
                total_fy += fy
                
                # Moment about centroid (positive counterclockwise)
                dx = force.point[0] - centroid[0]
                dy = force.point[1] - centroid[1]
                moment = fx * dy - fy * dx
                total_moment += moment
        
        # Calculate resultant properties
        resultant_magnitude = math.sqrt(total_fx**2 + total_fy**2)
        
        # Eccentricity and application point
        if abs(total_fx) > 1e-12 or abs(total_fy) > 1e-12:
            # Eccentricity magnitude
            eccentricity = abs(total_moment) / (resultant_magnitude + 1e-12)
            
            # Application point (simplified - along virtual face)
            if abs(total_fx) > 1e-12:
                app_x = centroid[0] + total_moment / total_fx
                app_y = centroid[1]
            else:
                app_x = centroid[0]
                app_y = centroid[1] + total_moment / (total_fy + 1e-12)
        else:
            eccentricity = float('inf')
            app_x, app_y = centroid
        
        return InterfaceResultant(
            block_a=block_a,
            block_b=block_b,
            centroid=centroid,
            virtual_face_length=virtual_length,
            tangent=virtual_tangent,
            normal=virtual_normal,
            resultant_fx=total_fx,
            resultant_fy=total_fy,
            resultant_magnitude=resultant_magnitude,
            moment_z=total_moment,
            eccentricity=eccentricity,
            application_point=(app_x, app_y)
        )
    
    def analyze_all_interfaces(self) -> List[InterfaceResultant]:
        """Analyze all interfaces in the structure"""
        interfaces = self.find_interfaces()
        results = []
        
        print(f"\nFound {len(interfaces)} interfaces to analyze:")
        
        for block_a, block_b in interfaces:
            print(f"  Analyzing interface {block_a} ? {block_b}")
            try:
                result = self.analyze_interface(block_a, block_b)
                results.append(result)
            except Exception as e:
                print(f"    Error: {e}")
        
        return results
    
    def save_results(self, results: List[InterfaceResultant], output_file: str) -> None:
        """Save results to CSV file"""
        print(f"\nSaving results to: {output_file}")
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow([
                'BlockA', 'BlockB', 'Centroid_X', 'Centroid_Y', 'Virtual_Face_Length',
                'Tangent_X', 'Tangent_Y', 'Normal_X', 'Normal_Y',
                'Resultant_Fx', 'Resultant_Fy', 'Resultant_Magnitude',
                'Moment_Z', 'Eccentricity', 'Application_Point_X', 'Application_Point_Y'
            ])
            
            # Data rows
            for r in results:
                writer.writerow([
                    r.block_a, r.block_b,
                    f"{r.centroid[0]:.6f}", f"{r.centroid[1]:.6f}",
                    f"{r.virtual_face_length:.6f}",
                    f"{r.tangent[0]:.6f}", f"{r.tangent[1]:.6f}",
                    f"{r.normal[0]:.6f}", f"{r.normal[1]:.6f}",
                    f"{r.resultant_fx:.3f}", f"{r.resultant_fy:.3f}",
                    f"{r.resultant_magnitude:.3f}",
                    f"{r.moment_z:.6f}", f"{r.eccentricity:.6f}",
                    f"{r.application_point[0]:.6f}", f"{r.application_point[1]:.6f}"
                ])
        
        print(f"Saved {len(results)} interface analyses")
    
    def print_summary(self, results: List[InterfaceResultant]) -> None:
        """Print a summary of all interfaces"""
        print("\n" + "="*80)
        print("INTERFACE RESULTANT SUMMARY")
        print("="*80)
        
        print(f"{'Interface':<12} {'Fx':<10} {'Fy':<10} {'|F|':<10} {'Moment':<12} {'Ecc.':<10}")
        print("-" * 80)
        
        for r in results:
            interface_name = f"{r.block_a}?{r.block_b}"
            print(f"{interface_name:<12} {r.resultant_fx:<10.1f} {r.resultant_fy:<10.1f} "
                  f"{r.resultant_magnitude:<10.1f} {r.moment_z:<12.3f} {r.eccentricity:<10.4f}")

def main():
    """Main function"""
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python interface_analyzer.py <results_file.txt>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = "interface_resultants.csv"
    
    try:
        # Create analyzer and parse file
        analyzer = InterfaceAnalyzer()
        analyzer.parse_results_file(input_file)
        
        # Analyze all interfaces
        results = analyzer.analyze_all_interfaces()
        
        # Print summary
        analyzer.print_summary(results)
        
        # Save to CSV
        analyzer.save_results(results, output_file)
        
        print(f"\nAnalysis complete! Results saved to {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()