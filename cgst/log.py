#!/usr/bin/env python3
"""
This module contains a class for writing output to both stdout and a log file.
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git
This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import textwrap


def log(message='', end='\n'):
    print(message, file=sys.stderr, flush=True, end=end)


def log_without_backline(message='', end=''):
    print(message, file=sys.stderr, flush=True, end=end)


def section_header(text):
    log()
    print(bold_yellow_underline(text), file=sys.stderr, flush=True)


def tool_output_log(text):
    log_without_backline()
    print(bold_cyan(text), file=sys.stderr, flush=True)


def tool_error_log(text):
    log_without_backline()
    print(bold_red_underline(text), file=sys.stderr, flush=True)


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
YELLOW = '\033[93m'
CYAN = '\033[96m'
RED = '\033[41m'
DIM = '\033[2m'


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def bold_red_underline(text):
    return RED + BOLD + UNDERLINE + text + END_FORMATTING


def bold_cyan(text):
    return CYAN + BOLD + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def explanation(text, indent_size=4):
    """
    This function writes explanatory text to the screen. It is wrapped to the terminal width for
    stdout but not wrapped for the log file.
    """
    text = ' ' * indent_size + text
    terminal_width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(dim(line))
    log()


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size


def write_log(log_file_path, log_message):
    try:
        with open(log_file_path, 'w') as log_file:
            log_file.write(log_message)
    except IOError as e:
        return e


def log_process_with_output_file(process, log_message, log_file_path):
    output_msg = error_msg = ""
    while True:
        output = ""
        error = ""
        if process.stdout is not None:
            output = process.stdout.readline().decode("utf-8").rstrip()
        if process.stderr is not None:
            error = process.stderr.readline().decode("utf-8").rstrip()
        if output == '' and process.poll() is not None:
            break
        if output:
            tool_output_log(output)
            output_msg += output
        if error:
            tool_error_log(error)
            error_msg += error

    # log
    log_message = f"{log_message}\n{output_msg}\n{error_msg}"
    write_log(log_file_path, log_message)


def log_process(process, string_to_search_in_log):
    output_list = []
    while True:
        output = ""
        error = ""
        if process.stdout is not None:
            output = process.stdout.readline().decode("utf-8").rstrip()
        if process.stderr is not None:
            error = process.stderr.readline().decode("utf-8").rstrip()
        if output == '' and process.poll() is not None:
            break
        if output:
            if string_to_search_in_log.lower() in output.lower():
                output_list.append(output.strip().split("\t")[0])
            tool_output_log(output)
        if error:
            tool_error_log(error)
    return output_list
