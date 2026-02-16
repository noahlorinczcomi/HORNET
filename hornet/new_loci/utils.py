import platform


def callPlink():
    sys = platform.system()
    if sys == "Linux":
        call_ = "./plinkdir/linux/plink"
    elif sys == "Darwin":
        call_ = "./plinkdir/mac/plink"
    elif sys == "Windows":
        call_ = "plinkdir/windows/plink.exe"
    return call_


def callDelete():
    sys = platform.system()
    if (sys == "Linux") | (sys == "Darwin"):
        call_ = "rm"
    elif sys == "Windows":
        call_ = "del"
    return call_


def flatten_list(_2d_list):  # flatten list of lists into a single list
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
