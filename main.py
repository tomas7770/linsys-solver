# System of linear equations solver, using Gauss–Jordan elimination

# A linear system is represented as a list of lists
# Each list contains the coefficients of a linear equation; the last number is the constant;
# all equations must have the same amount of coefficients
# e.g. [[1, 2, 3, 0], [1, -2, 1, 0]] represents {x+2y+3z=0, x-2y+z=0}, or the augmented coefficient matrix:
# [1  2 3 | 0]
# [1 -2 1 | 0]


def sum_rows(system, row_n, scalar, row_n_res):
    """
    Sums a multiple of a row of the coefficient matrix of the system with another row,
    storing the result in the latter row.

    :param system: a linear system
    :param int row_n: index of a row
    :param int or float scalar: multiply row of index row_n by this scalar before summing
    :param int row_n_res: index of another row
    :return: system
    """
    row = system[row_n]
    scaled_row = [scalar*x for x in row]
    row_res = system[row_n_res]
    new_row_res = [row_res[i]+scaled_row[i] for i in range(len(row_res))]
    system[row_n_res] = new_row_res
    return system


def scale_row(system, row_n, scalar):
    """
    Multiplies a row of the coefficient matrix of the system by a scalar.

    :param system: a linear system
    :param int row_n: index of a row
    :param int or float scalar: multiply row of index row_n by this scalar
    :return: system
    """
    row = system[row_n]
    scaled_row = [scalar*x for x in row]
    system[row_n] = scaled_row
    return system


def find_nonzero_element(row):
    """
    Finds first non-zero element of a row.

    :param row: a row of a linear system
    :return: (int) element index; len(row) if nonexistent
    """
    first_nonzero_elem = 0
    for i in range(len(row)):
        if row[i] != 0:
            break
        first_nonzero_elem += 1
    return first_nonzero_elem


def sort_rows(system):
    """
    Sorts the rows of the coefficient matrix of the system by index of first
    non-zero element.

    :param system: a linear system
    :return: system
    """
    system.sort(key=lambda row: find_nonzero_element(row))
    return system


def row_echelon(system):
    """
    Turns the coefficient matrix of the system into row echelon form.

    :param system: a linear system
    :return: system
    """
    row_amount = len(system)
    columns = len(system[0])
    for i in range(row_amount-1):
        sort_rows(system)
        row = system[i]
        column = find_nonzero_element(row)
        if column > columns-1:
            # Found a row of zeroes
            break
        elem = row[column]
        for j in range(i+1, row_amount):
            row_next = system[j]
            elem_next = row_next[column]
            if elem_next == 0:
                continue
            scalar = -elem_next/elem
            sum_rows(system, i, scalar, j)
    sort_rows(system)
    return system


def reduced_row_echelon(system):
    """
    Turns the coefficient matrix of the system into reduced row echelon form.

    :param system: a linear system
    :return: system
    """
    row_echelon(system)
    columns = len(system[0])
    for i in range(len(system)-1, -1, -1):
        row = system[i]
        column = find_nonzero_element(row)
        if column > columns-1:
            # Found a row of zeroes
            continue
        elem = row[column]
        scale_row(system, i, 1/elem)
        for j in range(i-1, -1, -1):
            row_prev = system[j]
            elem_prev = row_prev[column]
            if elem_prev == 0:
                continue
            # Optimization: first element of "row" is always 1, so no need to divide by that
            scalar = -elem_prev
            sum_rows(system, i, scalar, j)
    return system


def assert_system(system):
    """
    Errors if a linear system is not in a valid form.

    :param system: a linear system
    :returns: True
    """
    def assert_row(row):
        if not isinstance(row, list):
            raise TypeError("equation must be a list")
        if not all(map(lambda x: isinstance(x, (int, float)), row)):
            raise ValueError("equation must be a list of numbers")
    if not isinstance(system, list):
        raise TypeError("system must be a list")
    if not system:
        raise ValueError("system must have at least 1 equation")
    if not all(map(lambda l: isinstance(l, list), system)):
        raise ValueError("system must be a list of lists")
    map(assert_row, system)
    columns = len(system[0])
    if columns < 2:
        raise ValueError("equations must have at least 2 coefficients")
    if not all(map(lambda row: len(row) == columns, system)):
        raise ValueError("equations must have the same amount of coefficients")
    return True


def solve(system):
    """
    Solves a system of linear equations.

    :param system: a linear system
    :return: (tuple) solution set
    """
    assert_system(system)
    columns = len(system[0])
    solution = ["X"+str(i) for i in range(columns-1)]
    reduced_row_echelon(system)
    for row in reversed(system):
        first_nonzero_i = find_nonzero_element(row)
        if first_nonzero_i == columns-1:
            # 0 = non-zero element; impossible system
            return ()
        elif first_nonzero_i == columns:
            # Row of zeroes, skip
            continue
        var_value = ""
        for i in range(first_nonzero_i+1, columns):
            coefficient = row[i]
            if i == columns-1:
                # The constant
                var_value += "{0:+}".format(coefficient)
            elif coefficient != 0:
                # Coefficient of a variable
                var_value += "{0:+}*".format(-row[i])+solution[i]
        solution[first_nonzero_i] = var_value
    return tuple(solution)
