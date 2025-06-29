#!/usr/bin/env python3
"""
File Content Inspector - Let's see exactly what's in your results file
"""

def inspect_file(filepath):
    """Read and show the raw file content around the sections we need"""
    
    print("FILE CONTENT INSPECTOR")
    print("=" * 50)
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        print(f"File size: {len(content)} characters")
        print(f"Total lines: {len(content.split())}")
        
        # Look for our target sections
        forces_header = "=== FaceID; Pt; fn; ft ==="
        geometry_header = "--- Face ID, BlockJ, BlockK, Pt1, Pt2, Normal, Tangent ---"
        
        # Find forces section
        forces_pos = content.find(forces_header)
        if forces_pos != -1:
            print(f"\n✅ Found forces header at position {forces_pos}")
            
            # Show 500 characters after the header
            section_start = forces_pos + len(forces_header)
            section_content = content[section_start:section_start + 500]
            
            print("RAW CONTENT after forces header:")
            print("=" * 40)
            print(repr(section_content))  # repr shows \n, \r, etc.
            print("=" * 40)
            
            # Show line by line
            lines = section_content.split('\n')[:10]  # First 10 lines
            print(f"\nFirst 10 lines after forces header:")
            for i, line in enumerate(lines):
                print(f"Line {i}: '{line}'")
        else:
            print("❌ Forces header not found!")
        
        # Find geometry section
        geom_pos = content.find(geometry_header)
        if geom_pos != -1:
            print(f"\n✅ Found geometry header at position {geom_pos}")
            
            # Show 500 characters after the header
            section_start = geom_pos + len(geometry_header)
            section_content = content[section_start:section_start + 500]
            
            print("RAW CONTENT after geometry header:")
            print("=" * 40)
            print(repr(section_content))  # repr shows \n, \r, etc.
            print("=" * 40)
            
            # Show line by line
            lines = section_content.split('\n')[:10]  # First 10 lines
            print(f"\nFirst 10 lines after geometry header:")
            for i, line in enumerate(lines):
                print(f"Line {i}: '{line}'")
        else:
            print("❌ Geometry header not found!")
        
        # Show all section headers in the file
        print(f"\n\nALL SECTION HEADERS FOUND:")
        print("-" * 40)
        
        # Look for lines starting with === or ---
        lines = content.split('\n')
        for i, line in enumerate(lines):
            line = line.strip()
            if line.startswith('===') or line.startswith('---'):
                print(f"Line {i+1}: {line}")
                # Show next few lines
                for j in range(1, 4):
                    if i+j < len(lines):
                        next_line = lines[i+j].strip()
                        if next_line:
                            print(f"  +{j}: '{next_line}'")
                        else:
                            print(f"  +{j}: (empty)")
                print()
        
    except Exception as e:
        print(f"Error reading file: {e}")

def main():
    """Main function"""
    
    # Update this path to match your file
    input_file = r"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\results_cairo.txt"
    
    inspect_file(input_file)
    
    input("\nPress Enter to exit...")

if __name__ == "__main__":
    main()