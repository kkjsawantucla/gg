#Code generated using ChatGPT "lol"
import re
from collections import defaultdict

def parse_chemical_formula(formula):
    # Regular expression to match elements and their counts
    pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    
    # Dictionary to store elements and their counts
    elements = defaultdict(int)
    
    # Find all matches in the formula
    for match in pattern.finditer(formula):
        element, count = match.groups()
        # If no count is specified, it is 1
        count = int(count) if count else 1
        elements[element] += count
    
    return dict(elements)

def construct_matrix_A(ref):
    # Dictionary to hold the matrix A with elements as keys
    A = defaultdict(lambda: [0] * len(ref))
    
    # Iterate over the ref dictionary to populate matrix A
    for i, formula in enumerate(sorted(ref.keys())):
        parsed_formula = parse_chemical_formula(formula)
        for element, count in parsed_formula.items():
            A[element][i] = count
    
    # Convert defaultdict to a regular dictionary
    return dict(A)


def solve_chemical_equations(chemical_dict_A, chemical_dict_B):
    # Get the list of elements involved in both A and B
    elements_A = list(chemical_dict_A.keys())
    elements_B = list(chemical_dict_B.keys())
    
    # Ensure the elements are in the same order for both A and B
    elements = sorted(set(elements_A).union(elements_B))
    # Create matrix A and vector B based on the ordered elements
    A = np.array([chemical_dict_A.get(element, [0] * len(next(iter(chemical_dict_A.values())))) for element in elements])
    B = np.array([chemical_dict_B.get(element, 0) for element in elements]).reshape(len(elements),1)
    
    # Solve for x
    solution, residuals, rank, s = np.linalg.lstsq(A, B, rcond=None)
    
    return solution

def get_ref_coeff(ref,product_formula):

    # Create matrix A and vector B to solve Ax = B
    A = construct_matrix_A(ref)
    B = parse_chemical_formula(product_formula)

    # Use np.linalg.lstsq to solve the linear equation
    solution = solve_chemical_equations(A, B)

    ref_coeff={}
    for i, formula in enumerate(sorted(ref.keys())):
        ref_coeff[formula] = round(solution[i][0])
    
    return ref_coeff