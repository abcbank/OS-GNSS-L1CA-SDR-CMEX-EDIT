%Function performs cold start acquisition on the collected "data". It
%searches for GPS L1C signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 20 ms of raw IF signal from the front-end..
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and
% implementation of an open-source BDS-3 B1C/B2a SDR receiver.
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

acqResults = acquisition_cmex( ...
    settings, ...
    m.data, ...
    length(m.data), ...
    length(settings.acqSatelliteList), ...
    max(settings.acqSatelliteList) ...
);

% acqResults = cmex_acquisition( ...
%     settings, ...
%     m.data, ...
%     length(m.data), ...
%     length(settings.acqSatelliteList), ...
%     max(settings.acqSatelliteList) ...
% );
