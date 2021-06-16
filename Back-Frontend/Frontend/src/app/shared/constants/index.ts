import { IDashboardItem, ISidebarItem, IStudyType } from '../model';
import { IObjectType } from '../model/investigation.model';

export const SIDEBAR: ISidebarItem[] = [
  {
    id: 'Dashboard',
    icon: '',
    title: 'Dashboard',
    href: 'doe',
  },
  {
    id: 'About',
    icon: '',
    title: 'About',
    href: '',
  },
  {
    id: 'Settings',
    icon: '',
    title: 'Settings',
    href: '',
  },
  {
    id: 'Source Code',
    icon: '',
    title: 'Source Code',
    href: '',
  },
];

export const DASHBOARD_ITEMS: IDashboardItem[] = [
  {
    title: 'DoE',
    description:
      'Click here to Start Designing an Experiment. Once your Experiment has been designed, select a date for sample collection, and sit back and relax. We will do the rest for you.',
    href: '/doe',
  },
  {
    title: 'Results Retrieval',
    description:
      'Check the results of your completed experiments. You can download them in csv or other compatible formats or you can use our analytics online.',
    href: '',
  },
  {
    title: 'Module 3',
    description:
      'Amet minim mollit non deserunt ullamco est sit aliqua dolor do amet sint. ',
    href: '',
  },
  {
    title: 'Module 4',
    description:
      'Amet minim mollit non deserunt ullamco est sit aliqua dolor do amet sint. ',
    href: '',
  },
];

export const STUDY_TYPES: IStudyType[] = [
  {
    studyId: 1,
    studyName: 'Observational Study',
  },
  {
    studyId: 2,
    studyName: 'Interventional Study',
  },
];

export const OBJECT_TYPES: IObjectType[] = [
  {
    id: 'MATERIAL_OBJECT',
    name: 'Material Object',
  },
  {
    id: 'BIO_MATERIAL_OBJECT',
    name: 'Biomaterial Object',
  },
];

export const PUBLICATIONS = [
  { value: 'Publication 1', label: 'Publication 1' },
  { value: 'Publication 2', label: 'Publication 2' },
  { value: 'Publication 3', label: 'Publication 3' },
  { value: 'Publication 4', label: 'Publication 4' },
  { value: 'Publication 5', label: 'Publication 5' },
];

export const ISA_DOCUMENT_LICENSES = [
  { value: 'License 1', label: 'License 1' },
  { value: 'License 2', label: 'License 2' },
  { value: 'License 3', label: 'License 3' },
  { value: 'License 4', label: 'License 4' },
  { value: 'License 5', label: 'License 5' },
];

export const CONTACTS = [
  { value: 'Contact 1', label: 'Contact 1' },
  { value: 'Contact 2', label: 'Contact 2' },
  { value: 'Contact 3', label: 'Contact 3' },
  { value: 'Contact 4', label: 'Contact 4' },
  { value: 'Contact 5', label: 'Contact 5' },
];

export const DESIGN_TYPES = [
  { value: 'Design Type 1', label: 'Design Type 1' },
  { value: 'Design Type 2', label: 'Design Type 2' },
  { value: 'Design Type 3', label: 'Design Type 3' },
  { value: 'Design Type 4', label: 'Design Type 4' },
  { value: 'Design Type 5', label: 'Design Type 5' },
];

export const FACTOR_TYPES = [
  { value: 'Factor Type 1', label: 'Factor Type 1' },
  { value: 'Factor Type 2', label: 'Factor Type 2' },
  { value: 'Factor Type 3', label: 'Factor Type 3' },
  { value: 'Factor Type 4', label: 'Factor Type 4' },
  { value: 'Factor Type 5', label: 'Factor Type 5' },
];

export const PROTOCOLS = [
  { value: 'Protocol 1', label: 'Protocol 1' },
  { value: 'Protocol 2', label: 'Protocol 2' },
  { value: 'Protocol 3', label: 'Protocol 3' },
  { value: 'Protocol 4', label: 'Protocol 4' },
  { value: 'Protocol 5', label: 'Protocol 5' },
];

export const PERFORMERS = [
  { value: 'Performer 1', label: 'Performer 1' },
  { value: 'Performer 2', label: 'Performer 2' },
  { value: 'Performer 3', label: 'Performer 3' },
  { value: 'Performer 4', label: 'Performer 4' },
  { value: 'Performer 5', label: 'Performer 5' },
];

export const MEASUREMENT_TYPE_ANNOTATIONS = [
  { value: 'MT1', label: 'MT1' },
  { value: 'MT2', label: 'MT2' },
  { value: 'MT3', label: 'MT3' },
  { value: 'MT4', label: 'MT4' },
  { value: 'MT5', label: 'MT5' },
];

export const TECHNOLOGY_TYPE_ANNOTATIONS = [
  { value: 'TT1', label: 'TT1' },
  { value: 'TT2', label: 'TT2' },
  { value: 'TT3', label: 'TT3' },
  { value: 'TT4', label: 'TT4' },
  { value: 'TT5', label: 'TT5' },
];
