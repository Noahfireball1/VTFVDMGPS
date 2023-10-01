function printText(option)

programStyle = 'green';
textStyle = '*green';

switch option
    case 1
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading selected configuration file\n')
    case 2
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading selected configuration file\n')
    case 3
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading Diamond DA-40 aircraft properties...\n')
    case 4
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading Diamond DA-40 aerodynamic properties...\n')
    case 5
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading 3-blade propeller properties...\n')
    case 6
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading noise models for FVDM and selected IMUs...\n')
    case 7
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Starting Simulation\n\n')
    case 8
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle, 'Deep Integration of a Flight Vehicle Dynamic Model in a Vector Tracking Receiver\n')
        cprintf(programStyle,'\t\t\t')
        cprintf(textStyle,'Author: Noah Miller, nsm0014@auburn.edu\n')
        cprintf(programStyle,'\t\t\t')
        cprintf(textStyle,'Version: 1.3\n')
        cprintf(programStyle,'\t\t\t')
        cprintf(textStyle,'Last Updated: 9/14/2023\n\n')

        %% Please Select a Valid Configuration File to Get started
        cprintf(programStyle,'[VTFVDMGPS]\t')
        cprintf(textStyle,'Please select a valid configuration file to get started...\n\n')

    case 9
        cprintf(programStyle,'\n[VTFVDMGPS]\t')
        cprintf(textStyle, 'Generating satellite states for current date and time\n')

    case 10
        cprintf(programStyle,'\n[VTFVDMGPS]\t')
        cprintf(textStyle, 'Initializing simulations for Monte-Carlo analysis\n')
    case 11
        cprintf(programStyle,'\n[VTFVDMGPS]\t')
        cprintf(textStyle, 'Loading satellite states for selected configuration\n')
    case 12
        cprintf(programStyle,'\n[VTFVDMGPS]\t')
        cprintf(textStyle, 'No measurement update performed\n')

    otherwise
end

end