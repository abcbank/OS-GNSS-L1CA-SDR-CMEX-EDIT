function [trackResults, channel]= run_tracking_cmex(m, channel, settings)
% Performs GPS  L1C code and carrier tracking for all channels using
% data and pilot correlating approach.
%
%[trackResults, channel] = fullbandTracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.init
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

%CVS record:
%$Id: fullbandTracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

% Modified by Gyu-In Jee to track GPS L1C TMBOC(6,1,4/33) signal
% 2022/09/15

% -------- Number of acqusired signals ------------------------------------
TrackedNr =0 ;
for channelNr = 1:settings.numberOfChannels
    if channel(channelNr).status == 'T'
        TrackedNr = TrackedNr+1;
    end
end
%% Start processing channels ==============================================
[trackResults] = tracking_cmex(channel, m.Data, length(m.Data), settings);