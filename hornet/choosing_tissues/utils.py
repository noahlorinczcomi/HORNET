def flatten_list(_2d_list: list) -> list:
    """
    Docstring for flatten_list

    :param _2d_list: Description
    """
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list
