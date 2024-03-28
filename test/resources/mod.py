FILE_1 = "analytic_cube_solution_density1_old.txt"
FILE_1_OUT = "analytic_cube_solution_density1.txt"
FILE_42 = "analytic_cube_solution_density42_old.txt"
FILE_42_OUT = "analytic_cube_solution_density42.txt"

def mod(input_name, output_name):
    with open(input_name, 'r') as read_file, open(output_name, 'w') as write_file:
        lines = read_file.readlines()[1:]  # Skip the first line

        for line in lines:
            elements = line.split()  # Splits the line into multiple elements

            # Changes the sign of the last three elements
            elements[-1] = str(float(elements[-1]) * -1)
            elements[-2] = str(float(elements[-2]) * -1)
            elements[-3] = str(float(elements[-3]) * -1)

            write_file.write(' '.join(elements) + '\n')


if __name__ == '__main__':
    mod(FILE_42, FILE_42_OUT)
