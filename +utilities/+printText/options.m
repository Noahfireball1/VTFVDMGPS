function [outputArg1,outputArg2] = options(option)
programStyle = 'green';
textStyle = '*green';

switch option
    case 1
        utilities.cprintf(programStyle,'[VTFVDMGPS]\t')
        utilities.cprintf(textStyle, 'Rinex file for specified date not found. Downloading from https://cddis.nasa.gov\n')

    case 2
        utilities.cprintf(programStyle,'[VTFVDMGPS]\t')
        utilities.cprintf(textStyle, 'Rinex File successfully downloaded and unzipped!\n')

    case 3
        utilities.cprintf(programStyle,'[VTFVDMGPS]\t')
        utilities.cprintf(textStyle, 'https://cddis.nasa.gov username and password not found.\n\t\t\tPlease enter below:\n')

    otherwise
end
end

