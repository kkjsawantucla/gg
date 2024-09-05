# Code generated using ChatGPT "lol"
"""Import re"""

import re
from collections import defaultdict
import numpy as np
from ase.data import chemical_symbols

__author__ = "Kaustubh Sawant"


def parse_chemical_formula(formula):
    """
    Args:
        formula (str): chemical formula

    Returns:
        Dict: Break chemical formula into specific elements
    """
    # Regular expression to match elements and their counts
    pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
    # Dictionary to store elements and their counts
    elements = defaultdict(int)
    # Find all matches in the formula
    for match in pattern.finditer(formula):
        element, count = match.groups()
        # If no count is specified, it is 1
        count = int(count) if count else 1
        elements[element] += count
    return dict(elements)


def construct_matrix_a(ref):
    """
    Args:
        ref (Dict): Default dictionary with chemical formula and potential

    Returns:
        Dict: Dictionary broken down into specific elements
    """
    # Dictionary to hold the matrix A with elements as keys
    a = defaultdict(lambda: [0] * len(ref))
    # Iterate over the ref dictionary to populate matrix A
    for i, formula in enumerate(sorted(ref.keys())):
        parsed_formula = parse_chemical_formula(formula)
        for element, count in parsed_formula.items():
            a[element][i] = count
    # Convert default dict to a regular dictionary
    return dict(a)


def solve_chemical_equations(chemical_dict_a, chemical_dict_b):
    """ Solve Ax=B as dictionaries

    Args:
        chemical_dict_A (Dict):
        chemical_dict_B (Dict):

    Returns:
        np.array: Solution for Ax=B
    """
    # Get the list of elements involved in both A and B
    elements_a = list(chemical_dict_a.keys())
    elements_b = list(chemical_dict_b.keys())
    # Ensure the elements are in the same order for both A and B
    elements = sorted(set(elements_a).union(elements_b))
    # Create matrix A and vector B based on the ordered elements
    a = np.array(
        [
            chemical_dict_a.get(e, [0] * len(next(iter(chemical_dict_a.values()))))
            for e in elements
        ]
    )
    b = np.array([chemical_dict_b.get(e, 0) for e in elements]).reshape(
        len(elements), 1
    )
    # Solve for x
    solution, _, _, _ = np.linalg.lstsq(a, b, rcond=None)
    return solution


def get_ref_coeff(ref, product_formula):
    """
    Args:
        ref (Dict): reference chemical potential dictionary
        product_formula (str): chemical formula of the product

    Returns:
        Dict: coefficients for the reference chemical potential dictionary
    """

    # Create matrix A and vector B to solve Ax = B
    a = construct_matrix_a(ref)
    b = parse_chemical_formula(product_formula)

    # Use np.linalg.lstsq to solve the linear equation
    solution = solve_chemical_equations(a, b)

    ref_coeff = {}
    for i, formula in enumerate(sorted(ref.keys())):
        ref_coeff[formula] = solution[i][0]
    return ref_coeff


def is_element(symbol):
    """ Check if the symbol is an element
    Args:
        symbol (str):

    Returns:
        Bool:
    """
    return symbol in chemical_symbols


def is_chemical_formula(formula):
    """_summary_

    Args:
        formula (str):

    Returns:
        Bool:
    """
    # Regular expression to match chemical formulas
    pattern = r"^([A-Z][a-z]*\d*)+$"
    return bool(re.match(pattern, formula))

def parse_formula_to_list(formula):
    """ Convert chemical formula into a list
    Args:
        formula (str):

    Returns:
        list:
    """
    # Regular expression to match element symbols and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)

    # Create a list of individual symbols
    elements = []
    for (element, count) in matches:
        count = int(count) if count else 1
        elements.extend([element] * count)

    return elements
