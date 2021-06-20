import React from 'react';
import { ISidebarItem } from '../../../shared/model';
import ReactCollapsible from 'react-collapsible';

import './Sidebar.style.scss';
import { Link } from 'react-router-dom';

const SidebarItem = (props: {
  sidebarItem: ISidebarItem;
  itemClick: any;
  isActive: boolean;
}) => {
  const { sidebarItem, itemClick = () => {}, isActive } = props;
  const reactCollapsibleRef = React.createRef<ReactCollapsible>();
  const handleCollapsibleOpen = () => {
    (reactCollapsibleRef.current as any).innerRef.style.overflow = 'visible';
  };
  const handleCollapsibleClosing = () => {
    (reactCollapsibleRef.current as any).innerRef.style.overflow = 'hidden';
  };
  return (
    <div
      className={`sidebar__item ${
        isActive ? 'sidebar__item--active' : 'sidebar__item--inactive'
      }`}
      onClick={() => itemClick(sidebarItem)}
    >
      <ReactCollapsible
        trigger={
          <div className="sidebar__item__text collapsible-trigger--inner">
            {sidebarItem.title}
            {sidebarItem.children?.length
              ? [
                  <span className="sidebar-children-badge">
                    {sidebarItem.children?.length}
                  </span>,
                  <span className="material-icons trigger-icon">
                    keyboard_arrow_down
                  </span>,
                ]
              : null}
          </div>
        }
        onOpen={handleCollapsibleOpen}
        onClosing={handleCollapsibleClosing}
        ref={reactCollapsibleRef}
      >
        {sidebarItem.children?.length
          ? sidebarItem.children.map((sidebarItemChild, index) => (
              <SidebarItem
                key={index}
                isActive={false}
                sidebarItem={sidebarItemChild}
                itemClick={() => {}}
              />
            ))
          : null}
      </ReactCollapsible>
    </div>
  );
};

const Sidebar = (props: {
  sidebarItems?: ISidebarItem[];
  itemClick?: any;
  activeItem?: ISidebarItem;
  children?: React.ReactChildren | any;
}) => {
  const { sidebarItems = [], itemClick, activeItem, children } = props;
  return (
    <div className="sidebar--inner">
      <Link to="/">
        <div className="site-title">
          <span className="ellipse">B</span>
          <span className="site-title__text">BFAIR</span>
        </div>
      </Link>
      <div className="sidebar__items">
        {sidebarItems.map((sidebarItem: ISidebarItem, index) => (
          <SidebarItem
            key={index}
            isActive={activeItem?.id === sidebarItem.id}
            sidebarItem={sidebarItem}
            itemClick={itemClick}
          />
        ))}
      </div>
      {children}
    </div>
  );
};

export default Sidebar;
