import sys
import re

def find_line_number(file, include_line):
    """Find the line number of an include statement in a file."""
    with open(file, 'r') as f:
        for num, line in enumerate(f, 1):
            if include_line.strip() in line:
                return num
    return None

def main(iwyu_log):
    current_file = None
    include_pattern = re.compile(r'^\s*- #include ["<](.*)[">]')

    with open(iwyu_log, 'r') as log:
        file_lines = log.readlines()
    for line in file_lines:
        if line.endswith("should remove these lines:\n"):
            current_file = line.split()[0]
        elif include_pattern.match(line):
            include = include_pattern.match(line).groups()[0]
            line_number = find_line_number(current_file, include)
            if line_number:
                print(f"{line.strip()} // lines {line_number}-{line_number}")
                line_number = None
            else:
                print(f"WARNING: Could not find {include} in {current_file}", end='')
            continue
        print(line, end='')

if __name__ == "__main__":
    main(sys.argv[1])

