import React, { useState } from 'react';
import { connect } from 'react-redux';
import { ISidebarItem } from '../../../../../shared/model';
import { Investigation } from '../../../../../shared/model/investigation.model';
import Sidebar from '../../../../organisms/Sidebar';

const InvestigationFormSidebar = (props) => {
  const { investigation } = props;
  const InvestigationFormSidebarItems: ISidebarItem[] = [
    {
      id: investigation.id,
      icon: '',
      title: 'General Information',
      href: '#general',
    },
  ];
  const [activeItem, setActiveItem] = useState(
    InvestigationFormSidebarItems[0]
  );

  InvestigationFormSidebarItems.push(
    ...(investigation as Investigation).studies.map((study) => ({
      id: study.id,
      icon: '',
      title: study.title,
      href: '',
      children: study.assays.map((assay) => ({
        id: assay.id,
        title: assay.title,
        href: '',
        icon: '',
      })),
    }))
  );

  console.log(InvestigationFormSidebarItems);

  const itemClick = (item) => {
    console.log(item);
    setActiveItem(item);
  };
  return (
    <Sidebar
      sidebarItems={InvestigationFormSidebarItems}
      activeItem={activeItem}
      itemClick={itemClick}
    />
  );
};

const mapStateToProps = (state) => ({
  investigation: state,
});

export default connect(mapStateToProps)(InvestigationFormSidebar);
