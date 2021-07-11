import React from 'react';
import './Header.style.scss';

import notifications from '../../../../assets/images/bell_notifications.svg';
import profileImage from '../../../../assets/images/photo.png';

const Header = (props) => {
  const { pageTitle } = props;
  return (
    <header className="header">
      <div>
        <div className="page-title">{pageTitle}</div>
      </div>
      <div>
        <img alt="" src={notifications} />
        <div className="profile-detail">
          <img alt="" className="profile__image" src={profileImage} />
          <p className="profile__name">John Doe</p>
        </div>
      </div>
    </header>
  );
};

export default Header;
