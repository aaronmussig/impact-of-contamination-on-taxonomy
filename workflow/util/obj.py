def get_class_fqn(obj: object) -> str:
    """Return the fully qualified name of the class."""
    return f'{obj.__module__}.{obj.__class__.__name__}'


def is_int(obj) -> bool:
    try:
        int(obj)
        return True
    except ValueError:
        return False
