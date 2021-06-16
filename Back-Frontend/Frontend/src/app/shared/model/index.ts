export interface ISidebarItem {
  id: string;
  icon: string;
  title: string;
  href: string;
  children?: ISidebarItem[];
}

export interface IDashboardItem {
  href: string;
  title: string;
  description: string;
}

export interface IStudyType {
  studyName: string;
  studyId: string | number;
}
