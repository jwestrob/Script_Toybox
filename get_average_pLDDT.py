
import sys

def calculate_average_b_factor(file_path):
    total_b_factor = 0
    atom_count = 0
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Extracting the B-factor value from the second-to-last column
                try:
                    b_factor = float(line.split()[-2])
                except (ValueError, IndexError):
                    continue
                total_b_factor += b_factor
                atom_count += 1

    # Calculate the average B-factor
    if atom_count == 0:
        print(f"{file_path}: No ATOM records found.")
        return
    average_b_factor = total_b_factor / atom_count
    print(f"{file_path}: {average_b_factor:.2f}")

# Main execution
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python average_b_factor.py <path_to_pdb_file>")
    else:
        file_path = sys.argv[1]
        calculate_average_b_factor(file_path)
