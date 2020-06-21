# propkatrj --- https://github.com/Becksteinlab/propkatraj
# Copyright (c) 2013-2017 David Dotson, Ricky Sexton, Armin Zjajo, Oliver Beckstein
# Released under the GNU General Public License v3+

import os
from io import StringIO

import pandas as pd
import numpy as np

import propka.run as pk
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

import warnings
import logging


class ProkpaTraj(AnalysisBase):
    def __init__(self, atomgroup, ,skip_failure=False, **kwargs):
        """Init routine for the PropkaTraj class

        Parameters
        ----------
        atomgroup : :class:`MDAnalysis.Universe` or
                    :class:`MDAnalysis.AtomgGroup`
            Part of the system to obtain pKas for.
        skip_failure : bool
            If set to ``True``, skip frames where PROPKA fails. If ``False``
            raise an exception. The default is ``False``.
            Log file (at level warning) contains information on failed frames.

        Results
        -------
        pkas : :class:`pandas.DataFrame`
            DataFrame giving estimated pKa value for each residue for each
            trajectory frame. Residue numbers are given as column labels,
            times as row labels.

        """

        # store input atomgroup as atomgroup
        self.ag = atomgroup.atoms
        # temporary file to write to
        self.tmpfile = os.path.join(os.path.dirname(ag.universe.filename),
                                    'current.pdb')

    def _setup_frames(self, trajectory, start=None, stop=None, step=None):
        super(ProkpaTraj, self)._setup_frames(trajectory, start, stop, step)
        # Container to store pkas
        self.pkas = []
        self.num_failed_frames = 0
        self.failed_frames = []

    def _single_frame(self):
        pstream = mda.lib.util.NamedStream(StringIO(), self.tmpfile)
        self.ag.write(pstream)

        pstream.reset() # reset for reading

        try:
            mol = pk.single(pstream, optargs=['--quiet'])
        except (IndexError, AttributeError) as err:
            if not skip_failure:
                raise
            else:
                errmsg = f"failure on frame: {self._ts.frame}"
                warnings.warn(errmsg)
                self.num_failed_frames += 1
                self.failed_frames.append(self._ts.frame)
        finally:
            pstream.close(force=True)  # deallocate

        confname = mol.conformation_names[0]
        conformation = mol.conformations[confname]
        groups = conformation.get_titratable_groups()

        # extract pka estimates from each residue
        self.pkas.append([g.pka_value for g in groups])

    def _conclude(self):
        # Ouput failed frames
        if self.num_failed_frames > 0:
            perc_failure = (self.num_failed_frames / self.n_frames) * 100
            wmsg = (f"number of failed frames = {self.num_failed_frames}\n"
                    f"percentage failure = {perc_failure}")
                    f"failed frames: {self.failed_frames_log}")
            warnings.warn(wmsg)
        
        self.df = pd.DataFrame(pkas, index=pd.Float64Index(self.times, 
               name='time'), columns=[g.atom.resNumb for g in groups])

        
def get_propka(universe, sel='protein', start=None, stop=None, step=None, skip_failure=False):
    """Get and store pKas for titrateable residues near the binding site.
    Parameters
    ----------
    universe : :class:`MDAnalysis.Universe`
        Universe to obtain pKas for.
    sel : str, array_like
        Selection string to use for selecting atoms to use from given
        ``universe``. Can also be a numpy array or list of atom indices to use.
    start : int
        Frame of trajectory to start from. `None` means start from beginning.
    stop : int
        Frame of trajectory to end at. `None` means end at trajectory end.
    step : int
        Step by which to iterate through trajectory frames. propka is slow,
        so set according to how finely you need resulting timeseries.
    skip_failure : bool
        If set to ``True``, skip frames where PROPKA fails. If ``False``
        raise an exception. The default is ``False``.
        Log file (at level warning) contains information on failed frames.
    Results
    -------
    pkas : :class:`pandas.DataFrame`
        DataFrame giving estimated pKa value for each residue for each
        trajectory frame. Residue numbers are given as column labels, times as
        row labels.
    """
   
    # need AtomGroup to write out for propka
    if isinstance(sel, string_types):
        atomsel = universe.select_atoms(sel)
    elif isinstance(sel, (list, np.array)):
        atomsel = universe.atoms[sel]

    # "filename" for our stream
    # use same name so that propka overwrites
    newname = os.path.join(os.path.dirname(universe.filename), 'current.pdb')

    # progress logging output (because this is slow...)
    pm = mda.lib.log.ProgressMeter(universe.trajectory.n_frames,
                                   format="{step:5d}/{numsteps} t={time:12.3f} ps  "
                                   "[{percentage:5.1f}%]",
                                   interval=1)

    times = []
    pkas = []
    failed_frames = 0
    failed_frames_log = []
    for ts in universe.trajectory[start:stop:step]:
        pm.echo(ts.frame, time=ts.time)

        # we create a named stream to write the atoms of interest into
        pstream = mda.lib.util.NamedStream(StringIO(), newname)
        atomsel.write(pstream)

        pstream.reset()         # reset for reading

        # we feed the stream to propka, and it reads it as if it were a file on
        # disk
        try:
               mol = pk.single(pstream, optargs=['--quiet'])
        except (IndexError, AttributeError) as err:
                     #https://github.com/Becksteinlab/propkatraj/issues/13
                     #https://github.com/Becksteinlab/propkatraj/issues/10
               if not skip_failure:          
                   raise
               else:
                   err_msg = "{0} (failure {2}): failing frame {1}".format(
                         universe.trajectory.filename, ts.frame, failed_frames)
                   failed_frames += 1
                   failed_frames_log.append(ts.frame)
                   logging.warning(err_msg)
                   continue
        finally:
               pstream.close(force=True)  # deallocate

        # parse propka data structures to get out what we actually want
        confname = mol.conformation_names[0]
        conformation = mol.conformations[confname]
        groups = conformation.get_titratable_groups()

        # extract pka estimates from each residue
        pkas.append([g.pka_value for g in groups])

        # record time
        times.append(ts.time)

    if failed_frames_log:
       logging.warning('number of failed frames = {0}'.format(failed_frames))
       logging.warning('percent failure = {0:.3f}%'.format(float(failed_frames)/len(universe.trajectory)*100))
       logging.warning('failed frames: %r', failed_frames_log)
   
    # a `pandas.DataFrame` is a good data structure for this data
    df = pd.DataFrame(pkas, index=pd.Float64Index(times, name='time'),
                      columns=[g.atom.resNumb for g in groups])

    return df
