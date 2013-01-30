import os

if os.sys.platform == "win32":
    import msvcrt
    _poll = False
    
    def set_poll():
        global _poll
        _poll = True
    def set_normal():
        global _poll
        _poll = False
    def get_key():
        if _poll and not msvcrt.kbhit():
            return None
        return msvcrt.getch()
else:
    # from www.python.org/doc/faq/library, point 2.3
    import termios, fcntl
    original_terminal_attributes = termios.tcgetattr(os.sys.stdin.fileno())
    poll_terminal_attributes = termios.tcgetattr(os.sys.stdin.fileno())
    poll_terminal_attributes[0] &= ~termios.ICANON & ~termios.ECHO
    original_flags = fcntl.fcntl(os.sys.stdin.fileno(), fcntl.F_GETFL)
    def set_poll():
        termios.tcsetattr(os.sys.stdin.fileno(), termios.TCSANOW, poll_terminal_attributes)
        fcntl.fcntl(os.sys.stdin.fileno(), fcntl.F_SETFL, original_flags | os.O_NONBLOCK)
    def set_normal():
        termios.tcsetattr(os.sys.stdin.fileno(), termios.TCSANOW, original_terminal_attributes)
        fcntl.fcntl(os.sys.stdin.fileno(), fcntl.F_SETFL, original_flags)
    def get_key():
        try:
            return os.sys.stdin.read(1)
        except IOError:
            return None

