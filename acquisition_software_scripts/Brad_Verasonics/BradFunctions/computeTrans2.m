function Trans = computeTrans2(varargin)
% computeTrans - Assign array specifications for known transducer types
% Copyright 2001-2018 Verasonics, Inc.  All world-wide rights and remedies
% under all intellectual property laws and industrial property laws are
% reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% Arguments:
%    One argument:
%        1st arg - Trans
%        if Trans is a structure, return the full Trans structure. If ID string, return transducer name only.
%    Two arguments:
%        1st arg should be transducer name as string.
%        2nd arg should be parameter desired (eg 'maxHighVoltage').
%        Returns value of parameter.
%
% Use built-in defaults, but allow prior specification of some attributes, such as Trans.frequency.
%
%
% Trans structure field definitions: (Refer to the Vantage Sequence
% Programming Manual for more details)
%
% Trans.name: string value representing probe name as used by the
%    Vantage software.  Note name is case sensitive and must be spelled
%    exactly as listed in the "KnownTransducers" array!
% Trans.id: Matlab double set to a 24 bit unsigned integer value
%    representing the probe eeprom ID code as defined and interpreted by
%    Verasonics.  Note that the id value is often listed or displayed as a
%    string of six hex digits, but the Trans.id value must always be set to
%    the numeric equivalent of that hex string!
% Trans.units: string variable that must be set to either 'mm' or
%    'wavelengths'; set by default to 'mm' if not specified by user before
%    calling computeTrans.  The units selection applies to the following
%    Trans structure fields: elementWidth, ElementPos, lensCorrection.
% Trans.type: integer-value Matlab double representing one of the probe
%    geometry types supported by the Vantage software.  Four types are
%    currently defined:
%        type = 0: Linear or phased array; all elements in a line along the
%            X-axis.  Y and Z coordinates of all elements are set to zero
%        type = 1: Curved linear array; all elements are in the X Z plane.
%            Y coordinate of all elements is zero
%        type = 2: Three-dimensional array; each element can have an
%            arbitrary X Y Z position, with an arbitrary orientation.  All
%            five entries in Trans.ElementPos are used ( X Y Z position
%            coordinates and azimuth and elevation orientation angles).  A
%            2D array with all elements in the Z=0 plane and oriented
%            straight ahead is a special case of type 2.
%        type = 3: Annular array; each element is a full-circle annular
%            ring with all elements concentric around the Z-axis through
%            the center of the array.  Each element can have an arbitrary
%            width and diameter, orientation, and offset in the Z
%            dimension.
% Trans.connType: integer-value Matlab double representing an index to one
%    of the probe connector types supported by the system, and identifying
%    the specific connector type and pinout used by the probe.
% Trans.frequency: Matlab double set to the nominal center frequency of the
%    probe in MHz.  The system software uses this value to define the
%    wavelength that will used by all system parameters that are specified
%    in "wavelength" units.  Receive demodulation frequency can be
%    specified independently of Trans.frequency, but since the default is
%    to set Receive.demodF = Trans.frequency, specifying a hardware
%    supported frequency here is helpful.
% Trans.Bandwidth: A [1 X 2] Matlab double array, set to the approximate
%    lower and upper bandwidth limits of the probe in MHz at a cutoff level
%    of -3 dB (one way) or -6 dB (round trip).
% Trans.numelements: Integer value Matlab double set to the number of
%    individual elements in the probe.
% Trans.ElementPos: An N X 5 Matlab double array, where N is equal to
%    Trans.numelements.  For Trans.type = 0, 1, or 2 the first three
%    entries in each row are the X Y Z coordinates of the associated
%    element and the last two entries are the orientation angles in azimuth
%    and elevation.:
%        ElementPos row:  [ X  Y  Z  azimuth  elevation ]
%    If elevation orientation is always zero, it can be
%    omitted and an N X 4 array can be used.  For Trans.type = 3 (annular
%    arrays) the first two entries in each row are the radius to the inner
%    and outer edges of the element's annular ring; the 3rd and 4th entries
%    are the radius and Z coordinates of the equal-area center of the
%    element; and the 5th entry is the element's angular orientation with
%    respect to the Z axis.
%        ElementPos row:  [ ri ro rc zc  azimuth ]
%    All position coordinates are in either mm or
%    wavelength units as selected by Trans.units.  All angular orientation
%    entries are in radians.
% Trans.elementWidth: Nominal width of individual elements, in either mm or
%    wavelength units as selected by Trans.units.
% Trans.spacing and Trans.spacingMm: Nominal center-to-center spacing of
%    individual elements; spacingMm is always in mm units and spacing is
%    always in wavelength units regardless of the setting of Trans.units.
% Trans.lensCorrection: Effective one-way path length through the lens that
%    separates the face of the transducer from the active element, in
%    either mm or wavelength units as selected by Trans.units.
% Trans.ElementSens: A 1 X 101 array of Matlab doubles representing the
%    relative off-axis sensitivity of an element as a function of angle
%    from -pi/2 to pi/2 in steps of pi/100.
% Trans.elevationApertureMm: Matlab double representing the elevation
%    aperture in mm. Applies only to "1D" arrays (Trans.type 0 or 1).
% Trans.elevationFocusMm: Matlab double representing the elevation
%    focus depth in mm. Applies only to "1D" arrays (Trans.type 0 or 1).
% Trans.maxHighVoltage: Matlab double representing the maximum allowed
%    transmit voltage setting for the probe in Volts (peak, not peak-to-peak).
% Trans.impedance: N X 2 complex double array, where the first entry in
%    each row specifies a frequency in MHz and the second entry is the
%    complex impedance in Ohms of the probe elements as seen by the system
%    at that frequency. N can be any positive integer representing the
%    number of impedance vs frequency pairs that have been specified.
% Trans.radius and Trans.radiusMm: Optional fields that apply only to
%    Trans.type values of 1, 2, or 3.  If specified for type 1 they give
%    the radius of a circular arc in the X-Z plane on which all elements of
%    the curved linear array are located.  If specified for type 2 or 3
%    they give the radius of a spherical surface upon which all elements
%    of the array are located.  radiusMm is always in mm units and radius
%    is always in wavelength units regardless of the setting of
%    Trans.units.

% revised Apr 2018 for 3.4.2 and 3.5, to fix bugs and update parameters for
% the HIFUPlex probes, and for GE9LD.  See VTS-767, 769, 796.
% revised Jan 2018 to add HIFUPlex probes, new Trans.type = 3 for annular arrays,
%   and updates to Trans structure field definitions in the comments
% revised Dec 13, 2017 to add L22-14vX and L22-14vX-LF

% Known transducers and their corresponding ID, high voltage limit, and HVMux status. *** Do not modify these values. ***
KnownTransducers = {'L22-14v',    '02AB18',  30,   0};


% The probes listed above with the comment "no data" are only recognized in
% terms of the probe name and ID; there is no other Trans structure data
% provided by the computeTrans function for these probes.


switch nargin

    case 1
        Trans = varargin{1};
        if ~isstruct(Trans)  % if a structure is not provided as input, assume input is ID to translate into string.
            % input argument may be ID as a hex string, or the actual ID
            % numeric value
            if ischar(Trans)
                % convert hex string to number so it won't matter how many
                % leading zeros were provided
                probeID = hex2dec(Trans);
            else
                % not a string so presumably a numeric value
                probeID = Trans;
            end
            probeIDhex = num2str(probeID, '%06X');
            probenum = find(strcmpi(probeIDhex, KnownTransducers(:, 2)), 1);
            if isempty(probenum)
                Trans = 'Unknown';
            else
                Trans = KnownTransducers{probenum, 1}; % return the probe name (string value)
            end
            return
        end
        if ~isfield(Trans,'name'), error('computeTrans: Trans.name must be provided in input structure.'); end
        probenum = find(strcmpi(Trans.name, KnownTransducers(:, 1)), 1);
        if isempty(probenum), error('computeTrans: Trans.name not recognized as known transducer.'); end
        speedOfSound = 1.540;  % default speed of sound in mm/usec
        verbose = 2;
        if evalin('base','exist(''Resource'',''var'')&&isfield(Resource,''Parameters'')')
            if evalin('base','isfield(Resource.Parameters,''speedOfSound'')')
                speedOfSound = evalin('base','Resource.Parameters.speedOfSound')/1000; % speed of sound in mm/usec
            end
            if evalin('base','isfield(Resource.Parameters,''verbose'')')
                verbose = evalin('base','Resource.Parameters.verbose');
            else
                verbose = 2; % default value- display warnings and status messages
            end
        end

        % check for user-specified units, and print warning message if not found
        if ~isfield(Trans,'units') || isempty(Trans.units)
            if verbose > 0
                fprintf(2, 'Warning: Trans.units not specified; selecting default units of mm.\n');
                fprintf(2, 'If script requires wavelength units, add an explicit definition of\n');
                fprintf(2, '"Trans.units = ''wavelengths'';" before calling computeTrans.\n');
            end
            Trans.units = 'mm';
        end
        if ~strcmp(Trans.units, 'mm') && ~strcmp(Trans.units, 'wavelengths')
            error('computeTrans: Unrecognized value for Trans.units.  Must be ''mm'' or ''wavelengths''.');
        end

        % if Trans.frequency value has already been specified, we will use
        % it as is.  VSX and update() will confirm the value matches the
        % A/D sample rate constraints, and will exit with an error message
        % to the user if not.  Therefore we do not need to validate the
        % Trans.frequency value here (and could not, since we don't know
        % intended use of 4/3 sampling or interleave, etc.).
        if isfield(Trans,'frequency')
            if isempty(Trans.frequency)
                % if empty, remove it so cases below will assign default frequency
                Trans = rmfield(Trans, 'frequency');
            end
        end
        % also allow user-specified Bandwidth to override the default:
        if isfield(Trans,'Bandwidth')
            if isempty(Trans.Bandwidth)
                % if empty, remove it so cases below will assign default
                % Bandwidth
                Trans = rmfield(Trans, 'Bandwidth');
            end
        end

        Trans.lensCorrection = 0; % specify default value, in case it is not set for a particular transducer;
        Trans.id = hex2dec(KnownTransducers(probenum, 2));

        switch KnownTransducers{probenum, 1}
            case 'L22-14v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate),
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                    Trans.elevationApertureMm = 1.5; % active elevation aperture in mm
                    Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit
                
%             case 'L15-Xtech'
%                 if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
%                 % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate),
%                 % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
%                 % Manufacturer specified center frequency is 18.0 MHz +/- 10%
%                 if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [10, 18]; end
%                 Trans.type = 0; % Array geometry is linear (x values only)
%                 Trans.connType = 1; % Determined by the UTA being used################
%                 Trans.numelements = 128;
%                     Trans.elevationAperatureMm = 1.5; % active elevation focus depth from lens on face of transducer
%                     Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer
%                 Trans.elementWidth = 0.085; % element width in mm; assumes 20 micron kerf
%                 Trans.spacingMm = 0.100; % Spacing between elements in mm
%                 Trans.ElementPos = zeros(Trans.numelements,5);
%                 Trans.ElementPos(:,1) = Trans.spacingMm * (-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
%                 % Lens correction: from mfr data sheet, matching layers are 0.06 mm thick at 2145 m/sec
%                 % average velocity, and lens is 0.45 mm thick at 1147 m/sec
%                 Trans.lensCorrection = 1000 * speedOfSound * (0.06/2145 + 0.45/1147);
%                 % velocities in m/sec, result in mm
%                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit
           
            otherwise
                Trans.type = [];
                Trans.frequency = [];
                Trans.spacingMm = [];
                Trans.elementWidth = [];
                Trans.ElementSens = [];
                Trans.lensCorrection = [];
                Trans.ElementPos = [ 0 0 0 0 ];
                if verbose > 2
                    disp(' ');
                    disp(['computeTrans Status: Data not available for ', Trans.name]);
                    disp('Trans structure must be provided in user script.');
                    disp(' ');
                end
%                 fprintf(2, ['computeTrans: Data not available for ', Trans.name, ';\n']);
%                 fprintf(2, 'Trans structure must be provided in user script.\n');
%                 error(' ');

        end

        % Set a conservative value for maxHighVoltage, if not already defined
        if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

        % Now convert all units as required, based on Trans.units
        scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths

        % regardless of units, always provide spacing and radius in
        % wavelengths if they have been defined
        if isfield(Trans, 'spacingMm') && ~isempty(Trans.spacingMm)
            Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths.
        end
        if  isfield(Trans, 'radiusMm') && ~isempty(Trans.radiusMm)
            Trans.radius = Trans.radiusMm * scaleToWvl; % convert radiusMm to wavelengths
        end

        % define Trans.ElementSens based on Trans.elementWidth, but only if
        % user has not already defined it; assign default elementWidth if
        % it doesn't exist or is empty
        if ~isfield(Trans,'ElementSens')
            % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
            if ~isfield(Trans,'elementWidth') || isempty(Trans.elementWidth)
                % create default value of zero if not assigned in case
                % statements above (zero implies the element is a point
                % source)
                Trans.elementWidth = 0;
            end
            Theta = (-pi/2:pi/100:pi/2);
            Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
            % note at this point elementWidth is in mm, so we have to
            % convert to wavelengths for the ElementSens calculation
            eleWidthWl = Trans.elementWidth * scaleToWvl;
            if eleWidthWl < 0.01
                % avoid the divide by zero for very small values (in this
                % case the sinc function will be extremely close to 1.0 for
                % all Theta, so we only need the cos term)
                Trans.ElementSens = abs(cos(Theta));
            else
                Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
            end
        end


        if strcmp(Trans.units, 'wavelengths')
            % convert all mm unit variables to wavelengths.  Note columns 4
            % and 5 of the ElementPos array are angles in radians, and do
            % not require units conversion
            Trans.elementWidth = Trans.elementWidth * scaleToWvl;
            Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
            Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
            Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
            if Trans.type == 3
                % for type 3 annular arrays, the fourth column is also a distance, not an angle
                Trans.ElementPos(:,4) = Trans.ElementPos(:,4) * scaleToWvl;
            end
            Trans.lensCorrection = Trans.lensCorrection * scaleToWvl;
        end

    case 2
        % Two inputs provided - 1st input is Trans.name, 2nd is parameter desired (currently only
        % 'maxHighVoltage' allowed).
        nameStr = varargin{1};
        if ischar(nameStr)
            probenum = find(strcmpi(nameStr, KnownTransducers(:, 1)), 1);
            % if 1st input is not a recognized transducer name, set probenum to
            % zero to flag as unrecognized transducer
            if isempty(probenum)
                probenum = 0; % special case of zero will trigger assignment of default value
            end
            Param = varargin{2};
            switch Param
                case 'maxHighVoltage'
                    if probenum == 0
                        % unrecognized transducer name, so assign default
                        % maxHighVoltage limit at hw max value
                        Trans = 96;
                    else
                        % return maxHighVoltage from KnownTransducers array
                        Trans = KnownTransducers{probenum, 3};
                    end
                case 'HVMux'
                    if probenum == 0
                        % unrecognized transducer name, so assign default
                        % of no HVMux
                        Trans = 0;
                    else
                        % return HVMux status from KnownTransducers array
                        Trans = KnownTransducers{probenum, 4};
                    end
                otherwise
                    error('computeTrans: Unrecognized parameter in 2nd argument.');
            end
        else
            error('computeTrans: When called with 2 inputs, 1st input must be transducer name.');
        end

    otherwise
        error('computeTrans: computeTrans accepts one or two input arguments.');

end
return


