""" AVL Wrapper session and input classes
"""
import glob
import os
import shutil
import subprocess

from avlwrapper.config import default_config, IS_PYTHON_3
from avlwrapper.output import OutputReader

if IS_PYTHON_3:
    import tkinter as tk
    from tempfile import TemporaryDirectory

else:
    import Tkinter as tk
    import tempfile

    # simple class which provides TemporaryDirectory-esque functionality
    class TemporaryDirectory(object):
        def __init__(self, suffix='', prefix='', dir=None):
            self.name = tempfile.mkdtemp(suffix=suffix,
                                         prefix=prefix,
                                         dir=dir)

        def __enter__(self):
            return self.name

        def __exit__(self, exc, value, tb):
            self.cleanup()

        def cleanup(self):
            shutil.rmtree(self.name)


class Input(object):
    def to_string(self):
        raise NotImplementedError


class Parameter(Input):
    """Parameter used in the case definition"""
    def __init__(self, name, value, setting=None):
        """
        :param str name: Parameter name, if not in Case.CASE_PARAMETERS, it's
            assumed to by a control name
        :param float value: Parameter value
        :param str or None setting: Parameter setting,
            see `Case.VALID_SETTINGS`
        """
        self.name = name
        self.value = value

        # by default, a parameter is not constraint by as setting, but set to its own value
        if setting is None:
            self.setting = name
        else:
            self.setting = setting

    def to_string(self):
        return " {0:<12} -> {1:<12} = {2}\n".format(self.name,
                                                    self.setting,
                                                    self.value)


class State(Input):
    """State used in the case definition"""
    def __init__(self, name, value, unit=''):
        self.name = name
        self.value = value
        self.unit = unit

    def to_string(self):
        return " {0:<10} = {1:<10} {2}\n".format(self.name,
                                                 self.value,
                                                 self.unit)


class Case(Input):
    """AVL analysis case containing parameters and states"""
    CASE_PARAMETERS = {'alpha': 'alpha', 'beta': 'beta',
                       'roll_rate': 'pb/2V', 'pitch_rate': 'qc/2V',
                       'yaw_rate': 'rb/2V'}

    VALID_SETTINGS = {'alpha', 'beta', 'pb/2V', 'qc/2V',
                      'rb/2V', 'CL', 'CY', 'Cl', 'Cm', 'Cn',
                      'd1', 'd2', 'd3', 'd4', 'd5'}

    CASE_STATES = {'alpha': ('alpha', 0.0, 'deg'),
                   'beta': ('beta', 0.0, 'deg'),
                   'roll_rate': ('pb/2V', 0.0, ''),
                   'pitch_rate': ('qc/2V', 0.0, ''),
                   'yaw_rate': ('rb/2V', 0.0, ''),
                   'CL': ('CL', 0.0, ''),
                   'cd_p': ('CDo', None, ''),
                   'bank': ('bank', 0.0, 'deg'),
                   'elevation': ('elevation', 0.0, 'deg'),
                   'heading': ('heading', 0.0, 'deg'),
                   'mach': ('Mach', None, ''),
                   'velocity': ('velocity', 0.0, 'm/s'),
                   'density': ('density', 1.225, 'kg/m^3'),
                   'gravity': ('grav.acc.', 32.2, 'm/s^2'),
                   'turn_rad': ('turn_rad.', 0.0, 'm'),
                   'load_fac': ('load_fac.', 0.0, ''),
                   'X_cg': ('X_cg', None, 'm'),
                   'Y_cg': ('Y_cg', None, 'm'),
                   'Z_cg': ('Z_cg', None, 'm'),
                   'mass': ('mass', 1.0, 'kg'),
                   'Ixx': ('Ixx', 1.0, 'kg-m^2'),
                   'Iyy': ('Iyy', 1.0, 'kg-m^2'),
                   'Izz': ('Izz', 1.0, 'kg-m^2'),
                   'Ixy': ('Ixy', 0.0, 'kg-m^2'),
                   'Iyz': ('Iyz', 0.0, 'kg-m^2'),
                   'Izx': ('Izx', 0.0, 'kg-m^2'),
                   'visc_CL_a': ('visc CL_a', 0.0, ''),
                   'visc_CL_u': ('visc CL_u', 0.0, ''),
                   'visc_CM_a': ('visc CM_a', 0.0, ''),
                   'visc_CM_u': ('visc CM_u', 0.0, '')}

    def __init__(self, name, **kwargs):
        """
        :param str name: case name

        :param kwargs: key-value pairs
            keys should be Case.CASE_PARAMETERS, Case.CASE_STATES or a control.
            values should be a numeric value or a Parameter object
        """
        self.name = name
        self.number = 1
        self.parameters = self._set_default_parameters()
        self.states = self._set_default_states()
        # print('PRINTING STATES: \n')
        # print(self.states)
        self.controls = []
        for key, value in kwargs.items():
            # if a parameter object is given, add to the dict
            if isinstance(value, Parameter):
                self.parameters[key] = value
            elif isinstance(value, State):
                self.states[key].value = value
                print('Hello')
            else:
                # if the key is an existing case parameter, set the value
                if key in self.CASE_PARAMETERS:
                    param_str = self.CASE_PARAMETERS[key]
                    self.parameters[param_str].value = value
                elif key in self.CASE_STATES:
                    self.states[key].value = value
                # if an unknown key-value pair is given,
                # assume its a control and create a parameter
                else:
                    param_str = key
                    self.controls.append(key)
                    self.parameters[param_str] = Parameter(name=param_str,
                                                           value=value)

    def _set_default_parameters(self):
        # parameters default to 0.0
        return {name: Parameter(name=name, setting=name, value=0.0)
                for _, name in self.CASE_PARAMETERS.items()}

    def _set_default_states(self):
        return {key: State(name=value[0], value=value[1], unit=value[2])
                for key, value in self.CASE_STATES.items()}

    def _check(self):
        self._check_parameters()
        self._check_states()

    def _check_states(self):
        # print("Valid parameter settings are: {0}".format(self.CASE_STATES))
        for key in self.states.keys():
            if key not in self.CASE_STATES:
                raise InputError("Invalid state variable: {0}"
                                 .format(key),"Valid parameter settings are: {0}".format(self.CASE_STATES))

    def _check_parameters(self):
        for param in self.parameters.values():
            if (param.setting not in self.VALID_SETTINGS
                    and param.setting not in self.controls):
                raise InputError("Invalid setting on parameter: {0}."
                                 .format(param.name),"Valid parameter settings are: {0}".format(self.VALID_SETTINGS))

    def to_string(self):
        self._check()

        # case header
        case_str = (" " + "-"*45 + "\n Run case {0:<2}:  {1}\n\n"
                    .format(self.number, self.name))

        # write parameters
        for param in self.parameters.values():
            case_str += param.to_string()

        case_str += "\n"

        # write cases
        for state in self.states.values():
            case_str += state.to_string()

        return case_str


class Session(object):
    """Main class which handles AVL runs and input/output"""
    OUTPUTS = {'Totals': 'ft', 'SurfaceForces': 'fn',
               'StripForces': 'fs', 'ElementForces': 'fe',
               'StabilityDerivatives': 'st', 'BodyAxisDerivatives': 'sb',
               'HingeMoments': 'hm', 'StripShearMoments': 'vm'}

    def __init__(self, geometry, cases=None, name=None, config=default_config):
        """
        :param avlwrapper.Geometry geometry: AVL geometry
        :param typing.Sequence[Case] cases: Cases to include in input files
        :param str name: session name, defaults to geometry name
        :param avlwrapper.Configuration config: (optional) dictionary
            containing setting
        """

        self.config = config

        self.geometry = geometry
        self.cases = self._prepare_cases(cases)
        self.name = name or self.geometry.name

        self._results = None

    def _prepare_cases(self, cases):

        # guard for cases=None
        if cases is None:
            return []

        # If not set, make sure XYZref, Mach and CD0 default to geometry input
        geom_defaults = {'X_cg': self.geometry.point[0],
                         'Y_cg': self.geometry.point[1],
                         'Z_cg': self.geometry.point[2],
                         'mach': self.geometry.mach,
                         'cd_p': self.geometry.cd_p}

        for idx, case in enumerate(cases):
            case.number = idx + 1
            for key, val in geom_defaults.items():
                if case.states[key].value is None:
                    case.states[key].value = val
        return cases

    @property
    def model_file(self):
        return self.name + '.avl'

    @property
    def case_file(self):
        return self.name + '.case'

    @property
    def requested_output(self):
        requested_outputs = {k for k, v in self.config['output'].items()
                             if v.lower() == 'yes'}
        lc_outputs = {k.lower(): (k, v) for k, v in self.OUTPUTS.items()}

        outputs = {}
        for output in requested_outputs:
            if output not in lc_outputs:
                raise InputError("Invalid output: {}".format(output))
            name, ext = lc_outputs[output]
            outputs[name] = ext
        return outputs

    def _write_geometry(self, target_dir):
        model_path = os.path.join(target_dir, self.model_file)
        with open(model_path, 'w') as avl_file:
            avl_file.write(self.geometry.to_string())

    def _copy_airfoils(self, target_dir):
        airfoil_names = self.geometry.get_external_airfoil_names()
        current_dir = os.getcwd()
        for airfoil in airfoil_names:
            airfoil_path = os.path.join(current_dir, airfoil)
            shutil.copy(airfoil_path, target_dir)

    def _write_cases(self, target_dir):
        # AVL is limited to 25 cases
        if len(self.cases) > 25:
            raise InputError('Number of cases is larger than '
                             'the supported maximum of 25.')

        case_file_path = os.path.join(target_dir, self.case_file)

        with open(case_file_path, 'w') as case_file:
            for case in self.cases:
                case_file.write(case.to_string())

    def _write_analysis_files(self, target_dir):
        self._write_geometry(target_dir)
        self._copy_airfoils(target_dir)
        if self.cases is not None:
            self._write_cases(target_dir)

    def run_avl(self, cmds, pre_fn, post_fn):
        with TemporaryDirectory(prefix='avl_') as working_dir:
            pre_fn(working_dir)

            process = self._get_avl_process(working_dir)
            process.communicate(input=cmds.encode())
            process.wait()

            ret = post_fn(working_dir)
        return ret

    def _get_cases_run_cmds(self, cases):
        cmds = "oper\n"
        for case in cases:
            cmds += "{0}\nx\n".format(case.number)
            for _, ext in self.requested_output.items():
                out_file = self._get_output_filename(case, ext)
                cmds += "{cmd}\n{file}\n".format(cmd=ext,
                                                 file=out_file)
        return cmds

    @property
    def _load_files_cmds(self):
        cmds = "load {0}\n".format(self.model_file)
        if self.cases:
            cmds += "case {0}\n".format(self.case_file)
        return cmds

    @property
    def _run_all_cases_cmds(self):
        cmds = self._load_files_cmds
        if self.cases:
            cmds += self._get_cases_run_cmds(self.cases)
        else:
            cmds += "oper\n"
            cmds += "x\n"
        cmds += "\nquit\n"
        return cmds

    def run_all_cases(self):
        results = self.run_avl(cmds=self._run_all_cases_cmds,
                               pre_fn=self._write_analysis_files,
                               post_fn=self._read_results)
        return results

    def _get_avl_process(self, working_dir):
        # guard for avl not being present on the system.
        # this used to be check at config read, but this allows
        # dynamic setting of the configuration
        
        if 'avl_bin' not in self.config.settings:
            raise FileNotFoundError("AVL not found or not executable,"
                    " check the configuration file")

        stdin = subprocess.PIPE
        stdout = open(os.devnull, 'w') if not self.config[
            'show_stdout'] else None

        # Buffer size = 0 required for direct stdin/stdout access
        return subprocess.Popen(args=[self.config['avl_bin']],
                                stdin=stdin,
                                stdout=stdout,
                                bufsize=0,
                                cwd=working_dir)

    def _read_results(self, target_dir):
        results = dict()
        for case in self.cases:
            results[case.name] = dict()
            for output, ext in self.requested_output.items():
                file_name = self._get_output_filename(case, ext)
                file_path = os.path.join(target_dir, file_name)
                reader = OutputReader(file_path=file_path)
                results[case.name][output] = reader.get_content()
        return results

    def _get_output_filename(self, case, ext):
        out_file = "{base}-{case}.{ext}".format(base=self.name,
                                                case=case.number,
                                                ext=ext)
        return out_file

    def show_geometry(self):
        with TemporaryDirectory(prefix='avl_') as working_dir:
            self._write_geometry(working_dir)
            cmds = self._show_geometry_cmds
            avl = self._get_avl_process(working_dir)
            run_with_close_window(avl, cmds)

    def _get_plot(self, target_dir, plot_name):
        if 'gs_bin' not in self.config.settings:
            raise Exception("Ghostscript should be installed"
                            " and enabled in the configuration file")
        gs = self.config.settings['gs_bin']
        in_file = os.path.join(target_dir, 'plot.ps')
        out_file = os.path.join(os.getcwd(), plot_name + '.png')
        cmd = [gs, '-dBATCH', '-dNOPAUSE',
               '-sDEVICE=png16m', '-sOutputFile="{}"'.format(out_file), in_file]
        subprocess.call(cmd)
        if '%d' in out_file:
            return glob.glob(out_file.replace('%d', '*'))
        else:
            return [out_file]

    def save_geometry_plot(self):
        plot_name = self.name + '-geometry'
        cmds = self._hide_plot_cmds
        cmds += self._show_geometry_cmds
        cmds += "h\n\n\nquit\n"
        return self.run_avl(cmds=cmds,
                            pre_fn=self._write_geometry,
                            post_fn=lambda d: self._get_plot(d, plot_name))

    @property
    def _hide_plot_cmds(self):
        return "plop\ng\n\n"

    @property
    def _show_geometry_cmds(self):
        cmds = "load {0}\n".format(self.model_file)
        cmds += "oper\ng\n"
        return cmds

    def show_trefftz_plot(self, case_number):
        with TemporaryDirectory(prefix='avl_') as working_dir:
            self._write_analysis_files(working_dir)
            cmds = self._load_files_cmds
            cmds += self._show_trefftz_case_cmds(case_number)
            avl = self._get_avl_process(working_dir)
            run_with_close_window(avl, cmds)

    def save_trefftz_plots(self):
        cmds = self._hide_plot_cmds
        cmds += self._load_files_cmds
        if self.cases:
            for idx in range(1, len(self.cases) + 1):
                cmds += self._show_trefftz_case_cmds(idx)
                cmds += "h\n\n"
        else:
            cmds += "oper\nx\nt\nh\n\n"
        cmds += "\n\nquit\n"

        return self.run_avl(cmds=cmds,
                            pre_fn=self._write_analysis_files,
                            post_fn=lambda d: self._get_plot(d, self.name + '-trefftz-%d'))

    @staticmethod
    def _show_trefftz_case_cmds(case_number):
        cmds = "oper\n"
        cmds += "{}\nx\n".format(case_number)
        cmds += "t\n"
        return cmds

    def export_run_files(self, path=None):
        if path is None:
            path = os.path.join(os.getcwd(), self.name)
        if not os.path.exists(path):
            os.mkdir(path)
        self._write_analysis_files(path)
        print("Input files written to: {}".format(path))


class _CloseWindow(tk.Frame):
    def __init__(self, on_open=None, on_close=None, master=None):
        # On Python 2, tk.Frame is an old-style class
        tk.Frame.__init__(self, master)

        # Make sure window is on top
        master.call('wm', 'attributes', '.', '-topmost', '1')
        self.pack()
        self._on_open = on_open
        self._on_close = on_close
        self.close_button = self.create_button()

    def create_button(self):
        # add quit method to button press
        def on_close_wrapper():
            if self._on_close is not None:
                self._on_close()
            top = self.winfo_toplevel()
            top.destroy()
        close_button = tk.Button(self, text="Close",
                                 command=on_close_wrapper)
        close_button.pack()
        return close_button

    def mainloop(self, n=0):
        if self._on_open is not None:
            self._on_open()
        tk.Frame.mainloop(self, n)


class InputError(Exception):
    pass


def run_with_close_window(avl, cmds):
    quit_cmd = '\n\nquit\n'
    tk_root = tk.Tk()

    def open_fn(): avl.stdin.write(cmds.encode())

    def close_fn():
        avl.stdin.write(quit_cmd.encode())
        avl.wait()

    app = _CloseWindow(on_open=open_fn, on_close=close_fn, master=tk_root)
    app.mainloop()
