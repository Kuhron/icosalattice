# decorator to mark functions as only wanting to be called once when a module is initialized
# so that we're not redoing work like recreating the icosahedron vertices and their adjacency


def constant_maker(constant_name):
    def function_decorator(f):
        def wrapper(*args, calling_to_create_constant, **kwargs):
            if not calling_to_create_constant:
                raise Exception(f"do not call this function; use the constant {constant_name} instead")
            else:
                return f(*args, **kwargs)
        return wrapper
    return function_decorator
