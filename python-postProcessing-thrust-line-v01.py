#!/usr/bin/env python3
"""
Interface Resultant Analyzer for Masonry Structures
Computes resultant forces and eccentricities across block interfaces

USAGE:
1. Update the file paths in the main() function
2. Run: python interface_analyzer.py
3. Results saved to interface_resultants.csv

REQUIREMENTS:
- Input file must contain sections:
  === FaceID; Pt; fn; ft ===
  --- Face ID, BlockJ, BlockK, Pt1, Pt2, Normal, Tangent ---
"""

import re
import math
import csv
import sys
import os
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
class MegaBlock:
    """Mapping of virtual blocks to real mega blocks"""
    mega_id: int
    virtual_block_ids: List[int]

@dataclass
class VirtualInterfaceResultant:
    """Resultant force and properties for a virtual block interface"""
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

@dataclass
class MegaBlockInterfaceResultant:
    """Resultant force for mega block interface (combination of virtual interfaces)"""
    mega_block_a: int
    mega_block_b: int
    virtual_interfaces: List[Tuple[int, int]]  # List of (virtualBlockA, virtualBlockB) pairs
    centroid: Tuple[float, float]
    virtual_face_length: float
    tangent: Tuple[float, float]
    normal: Tuple[float, float]
    resultant_fx: float
    resultant_fy: float
    resultant_magnitude: float
    moment_z: float
    eccentricity: float
    application_point: Tuple[float, float]
    num_virtual_interfaces: int

class InterfaceAnalyzer:
    """Main analyzer class"""
    
    def __init__(self):
        self.face_forces: List[FaceForce] = []
        self.face_geometries: Dict[int, FaceGeometry] = {}
        self.mega_blocks: Dict[int, MegaBlock] = {}
        self.mega_interfaces: Dict[Tuple[int, int], List[int]] = {}
    
    def parse_results_file(self, filepath: str) -> None:
        """Parse the results file to extract face forces and geometry"""
        print(f"Parsing results file: {filepath}")
        
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Parse face forces section
        self._parse_face_forces(content)
        
        # Parse face geometry section  
        self._parse_face_geometry(content)
        
        # Parse mega blocks section (optional)
        self._parse_mega_blocks(content)
        
        # Parse mega block interfaces section (optional)
        self._parse_mega_interfaces(content)
        
        print(f"Loaded {len(self.face_forces)} force entries")
        print(f"Loaded {len(self.face_geometries)} face geometries")
        print(f"Loaded {len(self.mega_blocks)} mega blocks")
        print(f"Loaded {len(self.mega_interfaces)} mega block interfaces")
    
    def _parse_face_forces(self, content: str) -> None:
        """Parse the === FaceID; Pt; fn; ft === section"""
        # Find the forces section - stop at the next === section
        pattern = r'=== FaceID; Pt; fn; ft ===\s*\n(.*?)(?=\n===|\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            raise ValueError("Could not find face forces section")
        
        forces_text = match.group(1).strip()
        print(f"DEBUG: Full forces section captured {len(forces_text)} characters")
        
        lines = [line.strip() for line in forces_text.split('\n') if line.strip()]
        print(f"DEBUG: Found {len(lines)} non-empty force lines")
        
        for line_num, line in enumerate(lines):
            if line_num < 5 or line_num % 10 == 0:  # Show first 5 and every 10th line
                print(f"DEBUG: Processing force line {line_num}: '{line}'")
                
            # Parse: faceId; ptX,ptY,ptZ; fn; ft
            parts = line.split(';')
            if len(parts) != 4:
                print(f"DEBUG: Skipping line {line_num} - wrong number of parts: {len(parts)}")
                continue
                
            try:
                face_id = int(parts[0].strip())
                
                # Parse coordinates
                coords_str = parts[1].strip()
                coords = [float(x.strip()) for x in coords_str.split(',')]
                
                fn = float(parts[2].strip())
                ft = float(parts[3].strip())
                
                if line_num < 5:  # Only show details for first few
                    print(f"DEBUG: Added force - Face {face_id}, fn={fn}, ft={ft}")
                
                self.face_forces.append(FaceForce(
                    face_id=face_id,
                    point=(coords[0], coords[1], coords[2]),
                    fn=fn,
                    ft=ft
                ))
            except (ValueError, IndexError) as e:
                print(f"DEBUG: Error parsing force line {line_num}: {e}")
                continue
    
    def _parse_face_geometry(self, content: str) -> None:
        """Parse the --- Face ID, BlockJ, BlockK, Pt1, Pt2, Normal, Tangent --- section"""
        # Find the geometry section - capture everything until end of file
        pattern = r'--- Face ID, BlockJ, BlockK, Pt1, Pt2, Normal, Tangent ---\s*\n(.*?)(?=\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            raise ValueError("Could not find face geometry section")
        
        geometry_text = match.group(1).strip()
        print(f"\nDEBUG: Full geometry section captured {len(geometry_text)} characters")
        
        lines = [line.strip() for line in geometry_text.split('\n') if line.strip()]
        print(f"DEBUG: Found {len(lines)} non-empty geometry lines")
        
        for line_num, line in enumerate(lines):
            if line_num < 5 or line_num % 10 == 0:  # Show first 5 and every 10th line
                print(f"DEBUG: Processing geometry line {line_num}: '{line}'")
            
            # Parse: faceId;blockJ;blockK;x1,y1,z1; x2,y2,z2; nx,ny; tx,ty
            parts = line.split(';')
            if len(parts) != 7:
                print(f"DEBUG: Skipping line {line_num} - wrong number of parts: {len(parts)}")
                continue
            
            try:
                face_id = int(parts[0].strip())
                block_j = int(parts[1].strip())
                block_k = int(parts[2].strip())
                
                if line_num < 5:  # Only show details for first few
                    print(f"DEBUG: Face {face_id} connects blocks {block_j} ↔ {block_k}")
                
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
            except (ValueError, IndexError) as e:
                print(f"DEBUG: Error parsing geometry line {line_num}: {e}")
                continue
    
    def _parse_mega_blocks(self, content: str) -> None:
        """Parse the --- MegaBlock ID, Virtual Block IDs --- section"""
        # Find the mega blocks section (optional)
        pattern = r'--- MegaBlock ID, Virtual Block IDs ---\s*\n(.*?)(?=\n---|\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            print("DEBUG: No MegaBlocks section found (optional)")
            return
        
        mega_text = match.group(1).strip()
        print(f"\nDEBUG: Found MegaBlocks section with {len(mega_text)} characters")
        
        lines = [line.strip() for line in mega_text.split('\n') if line.strip()]
        print(f"DEBUG: Found {len(lines)} mega block definitions")
        
        for line_num, line in enumerate(lines):
            print(f"DEBUG: Processing mega block line {line_num}: '{line}'")
            
            # Parse: megaId,virtualId1;virtualId2;virtualId3
            parts = line.split(',')
            if len(parts) < 2:
                print(f"DEBUG: Skipping line {line_num} - insufficient parts")
                continue
            
            try:
                mega_id = int(parts[0].strip())
                
                # Parse virtual block IDs (handle both comma and semicolon separators)
                virtual_ids_str = parts[1].strip()
                
                # Skip if empty
                if not virtual_ids_str:
                    print(f"DEBUG: Skipping MegaBlock {mega_id} - empty virtual blocks field")
                    continue
                
                virtual_ids = []
                # Try both comma and semicolon separators
                if ';' in virtual_ids_str:
                    separators = virtual_ids_str.split(';')
                else:
                    separators = virtual_ids_str.split(',')
                
                for vid_str in separators:
                    vid_str = vid_str.strip()
                    if vid_str:
                        try:
                            vid = int(vid_str)
                            virtual_ids.append(vid)  # Accept ANY integer (including 0, -1, etc.)
                        except ValueError:
                            print(f"DEBUG: Skipping invalid virtual block ID '{vid_str}' in MegaBlock {mega_id}")
                
                if virtual_ids:  # Only add if we have valid virtual blocks
                    self.mega_blocks[mega_id] = MegaBlock(
                        mega_id=mega_id,
                        virtual_block_ids=virtual_ids
                    )
                    
                    print(f"DEBUG: Added MegaBlock {mega_id} with virtual blocks {virtual_ids}")
                else:
                    print(f"DEBUG: Skipping MegaBlock {mega_id} - no valid virtual block IDs")
                
            except (ValueError, IndexError) as e:
                print(f"DEBUG: Error parsing mega block line {line_num}: {e}")
                continue
    
    def _parse_mega_interfaces(self, content: str) -> None:
        """Parse the --- MegaBlockA, MegaBlockB, Face IDs --- section"""
        # Find the mega block interfaces section (optional)
        pattern = r'--- MegaBlockA, MegaBlockB, Face IDs ---\s*\n(.*?)(?=\n---|\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            print("DEBUG: No MegaBlockInterfaces section found (optional)")
            return
        
        interfaces_text = match.group(1).strip()
        print(f"\nDEBUG: Found MegaBlockInterfaces section with {len(interfaces_text)} characters")
        
        lines = [line.strip() for line in interfaces_text.split('\n') if line.strip()]
        print(f"DEBUG: Found {len(lines)} mega block interface definitions")
        
        for line_num, line in enumerate(lines):
            print(f"DEBUG: Processing mega interface line {line_num}: '{line}'")
            
            # Parse: megaBlockA,megaBlockB,faceId1,faceId2,faceId3,...
            parts = line.split(',')
            if len(parts) < 3:
                print(f"DEBUG: Skipping line {line_num} - insufficient parts (need at least megaA,megaB,faceId)")
                continue
            
            try:
                mega_a = int(parts[0].strip())
                mega_b = int(parts[1].strip())
                
                # Parse face IDs (all remaining parts)
                face_ids = []
                for i in range(2, len(parts)):
                    face_id_str = parts[i].strip()
                    if face_id_str:
                        face_ids.append(int(face_id_str))
                
                if face_ids:  # Only add if we have valid face IDs
                    # Always store with smaller mega block first for consistency
                    interface_key = (min(mega_a, mega_b), max(mega_a, mega_b))
                    self.mega_interfaces[interface_key] = face_ids
                    
                    print(f"DEBUG: Added mega interface {mega_a}↔{mega_b} with {len(face_ids)} faces: {face_ids}")
                else:
                    print(f"DEBUG: Skipping mega interface {mega_a}↔{mega_b} - no valid face IDs")
                
            except (ValueError, IndexError) as e:
                print(f"DEBUG: Error parsing mega interface line {line_num}: {e}")
                continue
    
    def find_interfaces(self) -> List[Tuple[int, int]]:
        """Find all unique block-to-block interfaces"""
        interfaces = set()
        
        print(f"\nDEBUG: Looking for interfaces in {len(self.face_geometries)} faces")
        
        for geom in self.face_geometries.values():
            print(f"DEBUG: Face {geom.face_id}: Block {geom.block_j} ↔ Block {geom.block_k}")
            
            # Only consider interfaces between structural blocks (positive IDs)
            if geom.block_j > 0 and geom.block_k > 0:
                # Always order as (smaller, larger) for consistency
                block_a = min(geom.block_j, geom.block_k)
                block_b = max(geom.block_j, geom.block_k)
                interfaces.add((block_a, block_b))
                print(f"DEBUG: Added interface {block_a} ↔ {block_b}")
            else:
                print(f"DEBUG: Skipping face {geom.face_id} - connects to support (block ≤ 0)")
        
        result = sorted(list(interfaces))
        print(f"DEBUG: Found {len(result)} unique interfaces: {result}")
        return result
    
    def get_virtual_to_mega_mapping(self) -> Dict[int, int]:
        """Create mapping from virtual block ID to mega block ID"""
        mapping = {}
        for mega_block in self.mega_blocks.values():
            for virtual_id in mega_block.virtual_block_ids:
                mapping[virtual_id] = mega_block.mega_id
        return mapping
    
    def find_mega_interfaces(self) -> List[Tuple[int, int]]:
        """Find all unique mega block interfaces"""
        if not self.mega_blocks:
            print("No mega blocks defined - skipping mega interface analysis")
            return []
        
        virtual_to_mega = self.get_virtual_to_mega_mapping()
        mega_interfaces = set()
        
        print(f"\nDEBUG: Virtual to Mega mapping: {virtual_to_mega}")
        
        # Check for unassigned virtual blocks
        all_virtual_blocks = set()
        for geom in self.face_geometries.values():
            if geom.block_j > 0:  # Only structural blocks
                all_virtual_blocks.add(geom.block_j)
            if geom.block_k > 0:
                all_virtual_blocks.add(geom.block_k)
        
        assigned_virtual_blocks = set(virtual_to_mega.keys())
        unassigned_blocks = all_virtual_blocks - assigned_virtual_blocks
        
        if unassigned_blocks:
            print(f"DEBUG: WARNING - These virtual blocks are not assigned to any mega block:")
            print(f"  Unassigned: {sorted(unassigned_blocks)}")
            print(f"  Assigned: {sorted(assigned_virtual_blocks)}")
            print(f"  Only interfaces between assigned blocks will create mega interfaces.")
        
        interfaces_skipped = 0
        interfaces_created = 0
        
        for geom in self.face_geometries.values():
            # Only consider interfaces between structural blocks (positive IDs)
            if geom.block_j > 0 and geom.block_k > 0:
                # Map virtual blocks to mega blocks
                mega_j = virtual_to_mega.get(geom.block_j)
                mega_k = virtual_to_mega.get(geom.block_k)
                
                if mega_j is not None and mega_k is not None and mega_j != mega_k:
                    # Interface between different mega blocks
                    mega_a = min(mega_j, mega_k)
                    mega_b = max(mega_j, mega_k)
                    mega_interfaces.add((mega_a, mega_b))
                    print(f"DEBUG: Face {geom.face_id} creates mega interface {mega_a} ↔ {mega_b}")
                    interfaces_created += 1
                elif mega_j is None or mega_k is None:
                    # One or both blocks not assigned to mega blocks
                    interfaces_skipped += 1
                # If mega_j == mega_k, it's internal to a mega block (no interface)
        
        result = sorted(list(mega_interfaces))
        print(f"DEBUG: Found {len(result)} unique mega interfaces: {result}")
        print(f"DEBUG: Skipped {interfaces_skipped} virtual interfaces (unassigned blocks)")
        print(f"DEBUG: Created {interfaces_created} mega interface connections")
        return result
    
    def analyze_mega_interface_direct(self, mega_a: int, mega_b: int, face_ids: List[int]) -> MegaBlockInterfaceResultant:
        """Analyze interface between two mega blocks using directly specified face IDs"""
        
        # Target mega block (forces act on the larger ID)
        target_mega = max(mega_a, mega_b)
        source_mega = min(mega_a, mega_b)
        
        print(f"DEBUG: Analyzing direct mega interface {source_mega}→{target_mega} with faces {face_ids}")
        
        # Get virtual blocks for target mega block
        if target_mega not in self.mega_blocks:
            raise ValueError(f"Target mega block {target_mega} not defined")
        target_virtual_blocks = set(self.mega_blocks[target_mega].virtual_block_ids)
        
        # Collect all vertices from interface faces for centroid calculation
        all_vertices = []
        valid_faces = []
        
        for face_id in face_ids:
            if face_id not in self.face_geometries:
                print(f"DEBUG: Warning - Face {face_id} not found in geometry, skipping")
                continue
                
            face = self.face_geometries[face_id]
            all_vertices.extend([face.pt1[:2], face.pt2[:2]])
            valid_faces.append(face)
        
        if not valid_faces:
            raise ValueError(f"No valid faces found for mega interface {mega_a}↔{mega_b}")
        
        # Calculate interface centroid
        unique_vertices = list(set(all_vertices))
        centroid_x = sum(v[0] for v in unique_vertices) / len(unique_vertices)
        centroid_y = sum(v[1] for v in unique_vertices) / len(unique_vertices)
        centroid = (centroid_x, centroid_y)
        
        # Find the longest face to determine interface orientation
        longest_face = max(valid_faces, key=lambda f: f.length)
        virtual_tangent = longest_face.edge_direction
        virtual_normal = (-virtual_tangent[1], virtual_tangent[0])
        virtual_length = longest_face.length * 2  # Multiply by 2 as requested
        
        print(f"DEBUG: Interface orientation from longest face {longest_face.face_id} (length={virtual_length:.3f})")
        
        # Sum forces with correct sign convention
        total_fx = 0.0
        total_fy = 0.0
        total_moment = 0.0
        
        for face in valid_faces:
            # Determine which virtual block from target mega block is involved
            target_virtual_block = None
            original_role = None
            
            if face.block_j in target_virtual_blocks:
                target_virtual_block = face.block_j
                original_role = "J"  # Was J in C# → need to flip sign
                sign = -1.0
            elif face.block_k in target_virtual_blocks:
                target_virtual_block = face.block_k  
                original_role = "K"  # Was K in C# → keep sign
                sign = +1.0
            else:
                print(f"DEBUG: Warning - Face {face.face_id} doesn't connect to target mega block {target_mega}")
                continue
            
            print(f"DEBUG: Face {face.face_id}: target_virtual={target_virtual_block}, role={original_role}, sign={sign}")
            
            # Get all forces for this face
            face_forces = [f for f in self.face_forces if f.face_id == face.face_id]
            
            for force in face_forces:
                # Apply forces with correct sign (replicating C# convention)
                fx = sign * (force.fn * face.normal[0] + force.ft * face.tangent[0])
                fy = sign * (force.fn * face.normal[1] + force.ft * face.tangent[1])
                
                total_fx += fx
                total_fy += fy
                
                # Moment about centroid (positive counterclockwise)
                dx = force.point[0] - centroid[0]
                dy = force.point[1] - centroid[1]
                moment = fx * dy - fy * dx
                total_moment += moment
        
        # Calculate application point correctly using moment equilibrium
        tangent_cross_force = total_fx * virtual_tangent[1] - total_fy * virtual_tangent[0]
        
        if abs(tangent_cross_force) > 1e-12:
            # Eccentricity along the interface tangent direction
            eccentricity_along_tangent = total_moment / tangent_cross_force
            
            # Application point = centroid + eccentricity * tangent_direction
            app_x = centroid[0] + eccentricity_along_tangent * virtual_tangent[0]
            app_y = centroid[1] + eccentricity_along_tangent * virtual_tangent[1]
        else:
            # Force is parallel to interface, no moment arm
            app_x, app_y = centroid
        
        # Calculate resultant magnitude
        resultant_magnitude = math.sqrt(total_fx**2 + total_fy**2)
        
        # Create list of virtual interfaces for compatibility
        virtual_interfaces = [(face.block_j, face.block_k) for face in valid_faces]
        
        return MegaBlockInterfaceResultant(
            mega_block_a=source_mega,
            mega_block_b=target_mega,
            virtual_interfaces=virtual_interfaces,
            centroid=centroid,
            virtual_face_length=virtual_length,
            tangent=virtual_tangent,
            normal=virtual_normal,
            resultant_fx=total_fx,
            resultant_fy=total_fy,
            resultant_magnitude=resultant_magnitude,
            moment_z=total_moment,
            eccentricity=0.0,
            application_point=(app_x, app_y),
            num_virtual_interfaces=len(valid_faces)
        )
    
    def analyze_mega_interface(self, mega_a: int, mega_b: int) -> MegaBlockInterfaceResultant:
        """Analyze interface between two mega blocks"""
        virtual_to_mega = self.get_virtual_to_mega_mapping()
        
        # Find all virtual interfaces that contribute to this mega interface
        virtual_interfaces = []
        interface_faces = []
        
        for geom in self.face_geometries.values():
            if geom.block_j > 0 and geom.block_k > 0:
                mega_j = virtual_to_mega.get(geom.block_j)
                mega_k = virtual_to_mega.get(geom.block_k)
                
                if ((mega_j == mega_a and mega_k == mega_b) or 
                    (mega_j == mega_b and mega_k == mega_a)):
                    virtual_interfaces.append((geom.block_j, geom.block_k))
                    interface_faces.append(geom)
        
        if not interface_faces:
            raise ValueError(f"No interface found between mega blocks {mega_a} and {mega_b}")
        
        print(f"DEBUG: Mega interface {mega_a}↔{mega_b} contains {len(virtual_interfaces)} virtual interfaces")
        
        # Collect all vertices from all interface faces
        all_vertices = []
        for face in interface_faces:
            all_vertices.extend([face.pt1[:2], face.pt2[:2]])
        
        # Remove duplicates and calculate centroid
        unique_vertices = list(set(all_vertices))
        centroid_x = sum(v[0] for v in unique_vertices) / len(unique_vertices)
        centroid_y = sum(v[1] for v in unique_vertices) / len(unique_vertices)
        centroid = (centroid_x, centroid_y)
        
        # Find the longest face to determine interface orientation
        longest_face = max(interface_faces, key=lambda f: f.length)
        virtual_tangent = longest_face.edge_direction
        virtual_normal = (-virtual_tangent[1], virtual_tangent[0])
        virtual_length = longest_face.length * 2  # Multiply by 2 as requested
        
        print(f"DEBUG: Interface orientation from longest face {longest_face.face_id} (length={virtual_length:.3f})")
        
        # Sum forces from all contributing faces
        total_fx = 0.0
        total_fy = 0.0
        total_moment = 0.0
        
        for face in interface_faces:
            # Determine which mega block this force acts on
            mega_j = virtual_to_mega.get(face.block_j)
            mega_k = virtual_to_mega.get(face.block_k)
            
            # Get all forces for this face
            face_forces = [f for f in self.face_forces if f.face_id == face.face_id]
            
            for force in face_forces:
                # Force acts on mega_b, determine sign
                if mega_k == mega_b:
                    sign = 1.0
                else:
                    sign = -1.0
                
                # Global force components
                fx = sign * (force.fn * face.normal[0] + force.ft * face.tangent[0])
                fy = sign * (force.fn * face.normal[1] + force.ft * face.tangent[1])
                
                total_fx += fx
                total_fy += fy
                
                # Moment about centroid
                dx = force.point[0] - centroid[0]
                dy = force.point[1] - centroid[1]
                moment = fx * dy - fy * dx
                total_moment += moment
        
        # Calculate application point correctly using moment equilibrium
        tangent_cross_force = total_fx * virtual_tangent[1] - total_fy * virtual_tangent[0]
        
        if abs(tangent_cross_force) > 1e-12:
            # Eccentricity along the interface tangent direction
            eccentricity_along_tangent = total_moment / tangent_cross_force
            
            # Application point = centroid + eccentricity * tangent_direction
            app_x = centroid[0] + eccentricity_along_tangent * virtual_tangent[0]
            app_y = centroid[1] + eccentricity_along_tangent * virtual_tangent[1]
        else:
            # Force is parallel to interface, no moment arm
            app_x, app_y = centroid
        
        # Calculate resultant magnitude
        resultant_magnitude = math.sqrt(total_fx**2 + total_fy**2)
        
        return MegaBlockInterfaceResultant(
            mega_block_a=mega_a,
            mega_block_b=mega_b,
            virtual_interfaces=virtual_interfaces,
            centroid=centroid,
            virtual_face_length=virtual_length,
            tangent=virtual_tangent,
            normal=virtual_normal,
            resultant_fx=total_fx,
            resultant_fy=total_fy,
            resultant_magnitude=resultant_magnitude,
            moment_z=total_moment,
            eccentricity=0.0,
            application_point=(app_x, app_y),
            num_virtual_interfaces=len(virtual_interfaces)
        )
    
    def analyze_interface(self, block_a: int, block_b: int) -> VirtualInterfaceResultant:
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
        
        # Find the longest face to determine interface orientation
        longest_face = max(interface_faces, key=lambda f: f.length)
        virtual_tangent = longest_face.edge_direction
        virtual_length = longest_face.length * 2  # Multiply by 2 as requested
        
        # Virtual face normal (perpendicular to tangent)
        virtual_normal = (-virtual_tangent[1], virtual_tangent[0])
        
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
        
        # Calculate application point correctly using moment equilibrium
        tangent_cross_force = total_fx * virtual_tangent[1] - total_fy * virtual_tangent[0]
        
        if abs(tangent_cross_force) > 1e-12:
            # Eccentricity along the interface tangent direction
            eccentricity_along_tangent = total_moment / tangent_cross_force
            
            # Application point = centroid + eccentricity * tangent_direction
            app_x = centroid[0] + eccentricity_along_tangent * virtual_tangent[0]
            app_y = centroid[1] + eccentricity_along_tangent * virtual_tangent[1]
            eccentricity = abs(eccentricity_along_tangent)
        else:
            # Force is parallel to interface, no moment arm
            app_x, app_y = centroid
            eccentricity = 0.0
        
        return VirtualInterfaceResultant(
            block_a=block_a,
            block_b=block_b,
            centroid=centroid,
            virtual_face_length=virtual_length,
            tangent=virtual_tangent,
            normal=virtual_normal,
            resultant_fx=total_fx,
            resultant_fy=total_fy,
            resultant_magnitude=math.sqrt(total_fx**2 + total_fy**2),
            moment_z=total_moment,
            eccentricity=eccentricity,
            application_point=(app_x, app_y)
        )
    
    def analyze_all_interfaces(self) -> Tuple[List[VirtualInterfaceResultant], List[MegaBlockInterfaceResultant]]:
        """Analyze all interfaces in the structure - both virtual and mega block levels"""
        # Analyze virtual block interfaces
        virtual_interfaces = self.find_interfaces()
        virtual_results = []
        
        print(f"\nFound {len(virtual_interfaces)} virtual block interfaces to analyze:")
        
        for block_a, block_b in virtual_interfaces:
            print(f"  Analyzing virtual interface {block_a} ↔ {block_b}")
            try:
                result = self.analyze_interface(block_a, block_b)
                virtual_results.append(result)
            except Exception as e:
                print(f"    Error: {e}")
        
        # Analyze direct mega block interfaces (new approach)
        mega_results = []
        
        if self.mega_interfaces:
            print(f"\nFound {len(self.mega_interfaces)} direct mega block interfaces to analyze:")
            
            for (mega_a, mega_b), face_ids in self.mega_interfaces.items():
                print(f"  Analyzing direct mega interface {mega_a} ↔ {mega_b} with faces {face_ids}")
                try:
                    result = self.analyze_mega_interface_direct(mega_a, mega_b, face_ids)
                    mega_results.append(result)
                except Exception as e:
                    print(f"    Error: {e}")
        else:
            # Fallback to old approach if no direct interfaces specified
            mega_interfaces = self.find_mega_interfaces()
            
            if mega_interfaces:
                print(f"\nFound {len(mega_interfaces)} inferred mega block interfaces to analyze:")
                
                for mega_a, mega_b in mega_interfaces:
                    print(f"  Analyzing inferred mega interface {mega_a} ↔ {mega_b}")
                    try:
                        result = self.analyze_mega_interface(mega_a, mega_b)
                        mega_results.append(result)
                    except Exception as e:
                        print(f"    Error: {e}")
            else:
                print("\nNo mega block interfaces found")
        
        return virtual_results, mega_results
    
    def save_results(self, virtual_results: List[VirtualInterfaceResultant], 
                     mega_results: List[MegaBlockInterfaceResultant], 
                     output_file: str) -> None:
        """Save mega interface results to CSV file"""
        
        # Only save mega interface results (skip virtual interfaces)
        if mega_results:
            print(f"\nSaving mega interface results to: {output_file}")
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)
                
                # Essential header for mega interfaces (removed tangent points)
                writer.writerow([
                    'MegaBlockA', 'MegaBlockB', 'Num_Virtual_Interfaces',
                    'Centroid_X', 'Centroid_Y', 'Virtual_Face_Length',
                    'Resultant_Fx', 'Resultant_Fy', 'Moment_Z', 
                    'Application_Point_X', 'Application_Point_Y'
                ])
                
                # Data rows for mega interfaces
                for r in mega_results:
                    writer.writerow([
                        r.mega_block_a, r.mega_block_b, r.num_virtual_interfaces,
                        f"{r.centroid[0]:.6f}", f"{r.centroid[1]:.6f}",
                        f"{r.virtual_face_length:.6f}",
                        f"{r.resultant_fx:.3f}", f"{r.resultant_fy:.3f}", 
                        f"{r.moment_z:.6f}",
                        f"{r.application_point[0]:.6f}", f"{r.application_point[1]:.6f}"
                    ])
            
            print(f"Saved {len(mega_results)} mega interface analyses")
        else:
            print("No mega interface results to save")
    
    def print_summary(self, virtual_results: List[VirtualInterfaceResultant], 
                     mega_results: List[MegaBlockInterfaceResultant]) -> None:
        """Print a summary focusing on mega block interfaces"""
        
        # Mega interfaces summary (primary focus)
        if mega_results:
            print("\n" + "="*80)
            print("MEGA BLOCK INTERFACE SUMMARY")
            print("="*80)
            
            print(f"{'Interface':<12} {'#Virt':<6} {'Fx':<10} {'Fy':<10} {'Moment':<12}")
            print("-" * 60)
            
            for r in mega_results:
                interface_name = f"M{r.mega_block_a}↔M{r.mega_block_b}"
                print(f"{interface_name:<12} {r.num_virtual_interfaces:<6} {r.resultant_fx:<10.1f} {r.resultant_fy:<10.1f} "
                      f"{r.moment_z:<12.3f}")
            
            # Show which virtual interfaces contribute to each mega interface
            print(f"\nMega Interface Composition:")
            for r in mega_results:
                virtual_list = ", ".join([f"{a}↔{b}" for a, b in r.virtual_interfaces])
                print(f"  M{r.mega_block_a}↔M{r.mega_block_b}: {virtual_list}")
        
        # Brief virtual interfaces summary (secondary)
        if virtual_results:
            print(f"\n📊 Virtual interface analysis: {len(virtual_results)} interfaces processed")
        
        if not virtual_results and not mega_results:
            print("\nNo interfaces found to analyze.")

def main():
    """Main function"""
    
    # ==========================================================================
    # CONFIGURATION: Edit these paths for your project
    # ==========================================================================
    input_file = r"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\results_cairo.txt"
    output_file = r"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\interface_resultants.csv"
    
    print("Interface Resultant Analyzer")
    print("="*50)
    print(f"Input file:  {input_file}")
    print(f"Output file: {output_file}")
    print("="*50)
    
    try:
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"❌ ERROR: Input file not found!")
            print(f"   Looking for: {input_file}")
            print(f"   Please check the path and make sure your C# solver has run.")
            input("Press Enter to exit...")
            return
        
        # Create analyzer and parse file
        analyzer = InterfaceAnalyzer()
        analyzer.parse_results_file(input_file)
        
        # Analyze all interfaces (both virtual and mega)
        virtual_results, mega_results = analyzer.analyze_all_interfaces()
        
        # Print summary
        analyzer.print_summary(virtual_results, mega_results)
        
        # Save to CSV files
        analyzer.save_results(virtual_results, mega_results, output_file)
        
        print(f"\n✅ Analysis complete!")
        print(f"📊 Summary: Analyzed {len(virtual_results)} virtual interfaces")
        if mega_results:
            print(f"📊 Summary: Analyzed {len(mega_results)} mega interfaces")
            print(f"   Results saved to: {output_file}")
        else:
            print(f"   No mega blocks defined - only virtual interface analysis completed")
            print(f"   (No output file created - mega interfaces only mode)")
        
        input("\nPress Enter to exit...")
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        print(f"\nThis might be caused by:")
        print(f"   - Missing or incorrect file format")
        print(f"   - Missing required sections in results file")
        print(f"   - File permission issues")
        input("\nPress Enter to exit...")
        sys.exit(1)

if __name__ == "__main__":
    main()